using Test

# Directory where the examples in Literate.jl format are
lectures_DIR = joinpath(splitpath(abspath(@__DIR__))[1:end-1], "src", "lectures")
lectures = [f for f in readdir(lectures_DIR) if startswith(f, "lecture")]

println("Testing source files used by Literate.jl")
@testset "$lecture" for lecture in lectures
    @testset "$f" for f in readdir(joinpath(lectures_DIR, lecture)) if endswith(f, ".jl")
        input = joinpath(lectures_DIR, lecture, f)
        include(f)
        @test true
    end
end
