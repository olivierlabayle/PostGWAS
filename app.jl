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

@app begin
    @in selected_locus = first(loci_ids)
    
    @out traces = []
    @out layout = PlotlyJS.Layout()

    @onchange selected_locus begin
      println(selected_locus)
      locus = loci_dataset[selected_locus]
      println("ok 0")
      positions = locus.POS
      log10pvalues = locus.LOG10P
      pips = locus.PIP
      variant_ids = locus.ID
      ld = locus.PHASED_R2
      cs = locus.CS
      println(unique(locus.CS))
      println("ok 1")

      subplots = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
      add_trace!(subplots,
          scatter(
          x=positions,
          y=log10pvalues,
          mode="markers",
          text=variant_ids,
          showlegend=false,
          marker = attr(
            color = ld,
            colorscale = "Viridis",
            showscale = true,
            colorbar = attr(
                title = "LD",
                y = 0.5,     # position higher up
                len = 0.5,   # roughly half the figure height
                yanchor = "bottom"
            ))
        ),
        row=1,
        col=1
      )
      add_trace!(subplots,
        scatter(
          x=positions,
          y=pips,
          showlegend=true,
          mode="markers",
          text=variant_ids,
          marker = attr(color = cs),
          name = "CS"
        ),
        row=2,
        col=1
      )
      relayout!(subplots; 
        xaxis2_title="Position",
        yaxis_title="-log₁₀(p-value)",
        yaxis2_title="PIP",
        show_legend=true,
        legend = attr(
          y = -0.1,  # centered below plot
          xanchor = "center",
          yanchor = "top"
      )
      )

      # shared_xaxes=true, vertical_spacing=0.02,
      # relayout!(subplots, showlegend=false, title_text=string("LocusZoom ",  selected_locus))
      println("ok 3")
      traces = subplots.plot.data
      layout = subplots.plot.layout
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
    cell([
      Stipple.select(
        :selected_locus,
        options = loci_ids,
        label = "Loci IDs",
        clearable = true,
        useinput = true,
      )
    ])
    cell([
      StipplePlotly.plot(:traces, layout=:layout, syncevents=true)
    ])
    ]
end

@page("/", ui)

end