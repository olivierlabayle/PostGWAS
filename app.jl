module App

using GenieFramework
using DataFrames
using PlotlyBase
using PlotlyJS
using StipplePlotly
using JLD2
@genietools

global loci_dataset = jldopen("/Users/olabayle/Dev/PostGWAS/rich_loci.jld2")
global loci_ids = collect(keys(loci_dataset))

function make_plot_component(locus)
	positions = collect(locus.POS)
	log10p = collect(locus.LOG10P)
	variant_ids = collect(locus.ID)
	ld = collect(locus.PHASED_R2)
	subplots = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
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
	relayout!(subplots;
		title = attr(
			text = string("Locus Zoom: ", first(locus.LOCUS_ID)),
			x = 0.5,
			xanchor = "center",
			font = attr(size=20, color="darkblue", family="Arial")
		),
		xaxis2_title="Position",
		yaxis_title="-log₁₀(p-value)",
		yaxis2_title="PIP",
		show_legend=true,
		legend=attr(
			title = attr(text = "CS"),
			y=0.5,
			yanchor="top"
		)
	)
	return subplots
end

function get_ensembl_consequences(locus, variant_id)
	variant_row = only(eachrow(locus[locus.ID .== variant_id, :]))
	consequences = variant_row[:FULL_ENSEMBL_ANNOTATIONS]["transcript_consequences"]
	# Not all dicts will have the same keys: we need to unify them to create a DataFrame from it
	all_columns = union((keys(d) for d in consequences)...)
	template_dict = Dict(key => missing for key in all_columns)
	return DataFrame([merge(template_dict, dict) for dict in consequences])
end

@app begin
	@in selected_locus_id = first(loci_ids)
	@in selected_variant = first(loci_ids)

	@out locus_variants = []
	@out locus_table = DataTable(DataFrame())
	@out ensembl_table = DataTable(DataFrame())

	@out traces = []
	@out layout = PlotlyJS.Layout()

	@onchange selected_locus_id begin
		selected_locus = loci_dataset[selected_locus_id]
		println(names(selected_locus))
		# update locus variants
		locus_variants = selected_locus.ID
		# update locus_table
		locus_table = DataTable(selected_locus[!, [:CHROM, :POS, :ID, :REF, :ALT, :ALT_FREQ, :MOST_SEVERE_CONSEQUENCE, :LOG10P, :BETA, :SE, :PIP, :CS]])
		println(selected_locus[1, "FULL_ENSEMBL_ANNOTATIONS"])
		# Update plots
		subplots = make_plot_component(selected_locus)
		traces = subplots.plot.data
		layout = subplots.plot.layout
	end

	@onchange selected_variant begin
		selected_locus = loci_dataset[selected_locus_id]
		get_ensembl_consequences(locus, selected_variant)
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
		cell([p("Welcome")])
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
		cell([GenieFramework.table(:locus_table, flat = true, bordered = true, title = "Variants Info")])

		# Locus Zoom Plot
		cell([
			StipplePlotly.plot(:traces, layout=:layout, syncevents=true)
		])

		# Variant Annotations
		cell([
			GenieFramework.select(:selected_variant, options = :locus_variants, label = "Select a variant")
		])
		cell([
			GenieFramework.table(:ensembl_table, flat = true, bordered = true, title = "ENSEMBL Consequences")
		])

	]
end

@page("/", ui)

end