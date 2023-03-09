using DataFrames, CSV, GLMakie, CairoMakie,Chain, DataFrameMacros, Statistics
GLMakie.activate!()
files = filter(contains(r"publication_model"),readdir())
df = CSV.read(files, DataFrame)

word_font = "Arial/Helvetica"
my_theme = Theme(
    Axis = (
        titlealign = :left,
        titlefont = word_font,
        xticklabelfont = word_font,
        yticklabelfont = word_font,
        ylabelfont = word_font,
        xlabelfont = word_font,
        xticklabelsize = 28,
        yticklabelsize = 28,
        xlabelsize = 30,
        ylabelsize = 30,
        titlesize = 25,
        titlegap = 10
    ),
    ColorBar =(labelsize = 20)
)
set_theme!(my_theme)

fig = Figure()
ga = fig[1:2, 1] = GridLayout()
gb = fig[1:2, 2] = GridLayout()
a = Axis(ga[2,1], ylabel = "Number of antagonists",xlabel = "Interspecific interaction multiplier")
b = Axis(gb[2,1],xlabel = "Interspecific interaction multiplier")
hideydecorations!(b)
a_data = @chain df begin
    @groupby(:alpha, :stressors)
    @combine(:value = mean(:bray))
end
b_data = @chain df begin
    @groupby(:alpha, :stressors)
    @combine(:value = mean(:kendall))
end
mn = minimum(vcat(a_data.value...,b_data.value...))
mx = maximum(vcat(a_data.value...,b_data.value...))
hm_1 = heatmap!(a,a_data.alpha,a_data.stressors,a_data.value)
hm_2 = heatmap!(b,b_data.alpha,b_data.stressors,b_data.value)
Colorbar(ga[1,1], hm_1,vertical = false,label = "Bray-Curtis Dissimilarity",
        labelsize = 20,labelfont = word_font )
Colorbar(gb[1,1], hm_2,vertical = false,label = "Kendall rank correlation" ,
            labelsize = 20, labelfont = word_font)


Label(ga[1,1, TopLeft()], "(a)",
    textsize = 36,
    font = word_font,
    padding = (0, 5, 5, 0),
    halign = :right)

Label(gb[1,1, TopLeft()], "(b)",
    textsize = 36,
    font = word_font,
    padding = (0, 5, 5, 0),
    halign = :right)

    t.plots[1][1][][1].fonts .= Makie.to_font(["Makie TeX Heros", "Makie TeX Heros Italic"])
notify(t.plots[1][1])

# Can tweak the plot manually, before switching to CairoMakie to save pdf

CairoMakie.activate!()
save("pub_fig2.png",fig)
