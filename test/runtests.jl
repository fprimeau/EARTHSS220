using Test

# Directory where the examples in Literate.jl format are
lectures_DIR = joinpath(splitpath(abspath(@__DIR__))[1:end-1]..., "src", "lectures")
lectures = [f for f in readdir(lectures_DIR) if startswith(f, "lecture")]

@testset "Testing source files used by Literate.jl" begin
    @testset "$lecture" for lecture in lectures
        source_files = [f for f in readdir(joinpath(lectures_DIR, lecture)) if endswith(f, ".jl")]
        @testset "$f" for f in source_files
            input = joinpath(lectures_DIR, lecture, f)
            include(input)
            @test true
        end
    end
end
