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

    @out positions = []
    @out log10pvalues = []
    
    @out traces = []
    @out plot_layout = PlotlyJS.Layout()

    @onchange selected_locus begin
      println(selected_locus)
      locus = loci_dataset[selected_locus]
      positions = locus.POS
      log10pvalues = locus.LOG10P
      traces = [scatter(
        x=positions,
        y=log10pvalues,
        mode="markers",
        text=locus.ID,
        name="Trace 1"
      )
      scatter(
        x=positions,
        y=locus.PIP,
        mode="markers",
        name="Trace 2"
      )
      ]
      plot_layout = PlotlyJS.Layout(
        title=string("Locus View: ", selected_locus),
        xaxis=attr(
            title="Position",
            showgrid=false
        ),
        yaxis=attr(
            title="-log10(p-value)",
            showgrid=true,
            range=[0, 20]
        )
      )
      println("ok")
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
      StipplePlotly.plot(:traces, layout=:plot_layout, syncevents=true)
    ])
    ]
end

@page("/", ui)

end