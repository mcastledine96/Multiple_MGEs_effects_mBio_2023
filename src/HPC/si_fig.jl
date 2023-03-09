cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DataFrames, CSV, CairoMakie,Chain, DataFrameMacros, Statistics
files = filter(contains(r"SI_[0-9]+_20.csv"),readdir())

word_font = "Arial/Helvetica"
my_theme = Theme(
    Axis = (
        titlealign = :left,
        titlefont = word_font,
        xticklabelfont = word_font,
        yticklabelfont = word_font,
        ylabelfont = word_font,
        xlabelfont = word_font,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xlabelsize = 20,
        ylabelsize = 20,
        titlesize = 20,
        titlegap = 20
    ),
    ColorBar =(labelsize = 20)
)
set_theme!(my_theme)


s = 0.1 + (parse(Int,ARGS[1]) ÷ 10)/10
    fs = []
for i in 1:90
    if 0.1 + (parse(Int,match(r"[0-9]+", files[i]).match) ÷ 10)/10 == s
        push!(fs,files[i])
    end
end

for (r,n) in enumerate(2:2:20)
    fig = Figure()
    fig2 = Figure()

    df = CSV.read(fs[1], DataFrame)
    df = df[df.n .== n,:]
    for i in 2:length(fs)
        df_i = CSV.read(fs[i], DataFrame)
        df_i = df_i[df_i.n .== n,:]
        df = vcat(df,df_i)
    end

    a = Axis(fig[2,1], ylabel = "Number of antagonists",xlabel = "Interspecific interaction multiplier")
    b = Axis(fig2[2,1],ylabel = "Number of antagonists",xlabel = "Interspecific interaction multiplier")
    df = df[.!isnan.(df.bray),:] 
    a_data = @chain df begin
        @groupby(:alpha, :stressors)
        @combine(:value = mean(:bray))
    end
    df = df[.!isnan.(df.kendall),:] 
    b_data = @chain df begin
        @groupby(:alpha, :stressors)
        @combine(:value = mean(:kendall))
    end

    hm_1 = heatmap!(a,a_data.alpha,a_data.stressors,a_data.value)
    hm_2 = heatmap!(b,b_data.alpha,b_data.stressors,b_data.value)

    Colorbar(fig[1,1], hm_1,vertical = false,label = "Bray-Curtis Dissimilarity: n = $n, σ = $s",
            labelsize = 20,labelfont = word_font )
    Colorbar(fig2[1,1], hm_2,vertical = false,label = "Kendall rank correlation: n = $n, σ = $s" ,
                labelsize = 20, labelfont = word_font)


    save("si_fig_$(s)_$(n).png",fig)
    save("si_fig2_$(s)_$(n).png",fig2)
end
