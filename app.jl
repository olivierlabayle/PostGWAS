module App

using GenieFramework
using DataFrames
using PlotlyBase
using PlotlyJS
using StipplePlotly
using JLD2
using HTTP
using JSON

@genietools

const ENSEMBL_SERVER = "https://rest.ensembl.org"

const GLOBAL_LOCK = ReentrantLock();

global loci_dataset = jldopen("/Users/olabayle/Dev/PostGWAS/rich_loci.jld2")
global loci_ids = collect(keys(loci_dataset))

global cache = Dict(
	"loci_dataset" => loci_dataset,
	"loci_ids" => loci_ids
)

function get_or_set!(cache, locus_id)
	lock(GLOBAL_LOCK) do
		if haskey(cache, locus_id)
			return cache[locus_id]
		else
			cache[locus_id] = cache["loci_dataset"][locus_id]
		end
	end
end

function make_plot_component(locus)
	positions = collect(locus.POS)
	log10p = collect(locus.LOG10P)
	variant_ids = collect(locus.ID)
	ld = collect(locus.PHASED_R2)
	region = string(first(locus.CHROM), ":", minimum(locus.POS), "-", maximum(locus.POS))
    genes = get_genomic_features(region; features=["gene"])
	subplots = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.02)
	# GWAS Plot
	add_trace!(subplots,
		scatter(
			x=positions,
			y=log10p,
			mode="markers",
			text=variant_ids,
			showlegend=false,
			marker=attr(
				color=ld,
				colorscale="Viridis",
				showscale=true,
				colorbar=attr(
					title="LD",
					y=0.5,     # position higher up
					len=0.5,   # roughly half the figure height
					yanchor="bottom"
				)
			)
		),
		row=1,
		col=1
	)
	# Finemapping plot
	for (cs_key, cs_group) in pairs(groupby(locus[!, [:POS, :CS, :PIP, :ID]], :CS))
		group_positions = collect(cs_group.POS)
		group_pips = collect(cs_group.PIP)
		group_variant_ids = collect(cs_group.ID)
		marker_color = cs_key.CS isa Missing ? "grey" : cs_key.CS
		show_leg = cs_key.CS isa Missing ? false : true
		add_trace!(subplots,
			scatter(
				x=group_positions,
				y=group_pips,
				showlegend=show_leg,
				mode="markers",
				text=group_variant_ids,
				marker_color=marker_color,
				name=string(nrow(cs_group))
			),
			row=2,
			col=1
		)
	end
	# Plot Genes
	for (y_coord, feature) in enumerate(genes)
		x_range = feature["start"]:feature["end"]
		add_trace!(subplots,
			scatter(
				x=x_range,
				y=fill(y_coord, length(x_range)),
				mode="line",
				showlegend=false,
				line = attr(color = "blue", width = 5)
			),
			row=3,
			col=1
		)
		add_trace!(subplots,
			scatter(
				x=[(x_range[end]+x_range[1])/2],
				y=[y_coord+0.4],
				textposition = "top center",
				mode="text",
				text=[feature["gene_id"]],
				showlegend=false,
			),
			row=3,
			col=1
		)
    end
	# Layout update
	relayout!(subplots;
		title = attr(
			text = string("Locus Zoom: ", first(locus.LOCUS_ID)),
			x = 0.5,
			xanchor = "center",
			font = attr(size=20, color="darkblue", family="Arial")
		),
		xaxis3_title="Position",
		yaxis_title="-log₁₀(p-value)",
		yaxis2_title="PIP",
		yaxis3_title="Genes",
		show_legend=true,
		legend=attr(
			title = attr(text = "CS"),
			y=0.5,
			yanchor="top"
		)
	)
	return subplots
end

function get_ensembl_consequences(ensembl_annotations)
	consequences = ensembl_annotations["transcript_consequences"]
	# Not all dicts will have the same keys: we need to unify them to create a DataFrame from it
	all_columns = union((keys(d) for d in consequences)...)
	template_dict = Dict(key => missing for key in all_columns)
	return DataTable(DataFrame([merge(template_dict, dict) for dict in consequences]))
end

function get_genomic_features(region; features=["gene", "transcript", "cds", "exon", "regulatory", "motif"])
    ext = string("/overlap/region/human/", region, "?", join(map(f -> "feature=$(f)", features), ";"))
    headers=Dict("Content-Type" => "application/json", "Accept" => "application/json")
    r = HTTP.get(ENSEMBL_SERVER*ext, headers)
    return JSON.parse(String(r.body))
end

get_gtex_annotations(info::Missing) = DataTable(DataFrame())

get_gtex_annotations(info) = DataTable(DataFrame(info))

make_locus_table(locus) = 
	DataTable(locus[!, [:CHROM, :POS, :ID, :REF, :ALT, :ALT_FREQ, :MOST_SEVERE_CONSEQUENCE, :LOG10P, :BETA, :SE, :PIP, :CS]])

@app begin
	@in selected_locus_id = ""
	@in selected_variant = ""

	@out locus_variants = []
	@out locus_table = DataTable(DataFrame())

	@out ensembl_table = DataTable(DataFrame())
	@in ensembl_table_search = ""

	@out eqtl_table = DataTable(DataFrame())
	@in eqtl_table_search = ""

	@out sqtl_table = DataTable(DataFrame())
	@in sqtl_table_search = ""

	@out traces = []
	@out layout = PlotlyJS.Layout()

	@onchange selected_locus_id begin
		@info string("Updating locus: ", selected_locus_id)
		selected_locus = get_or_set!(cache, selected_locus_id)
		# update locus variants
		locus_variants = selected_locus.ID
		# update locus_table
		locus_table = make_locus_table(selected_locus)
		# Update plots
		@info string("Before plot")
		subplots = make_plot_component(selected_locus)
		@info string("After plot")
		println(length(subplots.plot.data))
		traces = subplots.plot.data
		layout = subplots.plot.layout
		# Update selected variant
		selected_variant = selected_locus_id
	end

	@onchange selected_variant begin
		@info string("Updating selected variant: ", selected_variant)
		selected_locus = get_or_set!(cache, selected_locus_id)
		variant_row = only(eachrow(selected_locus[selected_locus.ID .== selected_variant, :]))
		ensembl_table = get_ensembl_consequences(variant_row.FULL_ENSEMBL_ANNOTATIONS)
		eqtl_table = get_gtex_annotations(variant_row.GTEX_EQTL_INFO)
		sqtl_table = get_gtex_annotations(variant_row.GTEX_SQTL_INFO)
	end

	# Reactive plots interactions
	@mixin traces::PlotlyEvents
	@onchange traces_click begin
		@show traces_click
	end
	@onchange traces_hover begin
		@show traces_hover
	end
	@onchange traces_selected begin
		@show traces_selected
	end
	@onchange traces_relayout begin
		@show traces_relayout
	end
end

@mounted watchplots()

function ui()
	[
		cell([
			header("LocusZoom",
				class="text-center", 
				style="font-size:24px; font-weight:bold;"
			); 
			p("Welcome to LocusZoom, this page provides an integrated view of GWAS, fine-mapping and annotation results for a given locus.");
			p("Here is what you can do:");
			ol([
				li("Select a locus from the drop-down list below."),
				li("Browse the locus' summary statistics."),
				li("Navigate GWAS and fine-mapping tracks."),
				li("Select a variant for which you would like to see more information")
			])
		])
		# Locus Selection
		cell([
			GenieFramework.select(
				:selected_locus_id,
				options=loci_ids,
				label="Select a locus ID",
				clearable=true,
				useinput=true,
			)
		])
		# Variants Info
		cell([
			header("Results for locus: {{selected_locus_id}}", 
				class="text-center", 
				style="font-size:24px; font-weight:bold;"
			),
			GenieFramework.table(:locus_table, 
				flat = true, 
				bordered = true, 
				title = "Summary Statistics"
			)
		])
		# Locus Zoom Plot
		cell([
			StipplePlotly.plot(:traces, layout=:layout, syncevents=true)
		])

		# Variant Annotations
		cell([header("Browse Variant Annotations", 
				class="text-center", 
				style="font-size:24px; font-weight:bold;")])
		cell([
			GenieFramework.textfield("Variant ID", :selected_variant)
		])
		row([
			cell(class="st-col col-4 col-sm st-module", [
				GenieFramework.table(:ensembl_table, 
					flat = true, 
					bordered = true, 
					title = "ENSEMBL Consequences",
					var"row-key" = "name",
					filter = :ensembl_table_search,
					template(
						var"v-slot:top-right" = "",
						textfield(
							"",
							:ensembl_table_search,
							dense = true,
							debounce = "300",
							placeholder = "Search",
							[template(var"v-slot:append" = true, icon("search"))],
						),
					)
			)
			]),
			cell(class="st-col col-4 col-sm st-module", [
				GenieFramework.table(:eqtl_table, 
					flat = true, 
					bordered = true, 
					title = "GTEx eQTLs",
					var"row-key" = "name",
					filter = :eqtl_table_search,
					template(
						var"v-slot:top-right" = "",
						textfield(
							"",
							:eqtl_table_search,
							dense = true,
							debounce = "300",
							placeholder = "Search",
							[template(var"v-slot:append" = true, icon("search"))],
						),
					)
				)
			]),
			cell(class="st-col col-4 col-sm st-module", [
				GenieFramework.table(:sqtl_table, 
					flat = true, 
					bordered = true, 
					title = "GTEx sQTLs",
					var"row-key" = "name",
					filter = :sqtl_table_search,
					template(
						var"v-slot:top-right" = "",
						textfield(
							"",
							:sqtl_table_search,
							dense = true,
							debounce = "300",
							placeholder = "Search",
							[template(var"v-slot:append" = true, icon("search"))],
						),
					)
				)
			])
		])

	]
end

@page("/", ui)

end