# generate examples
using Literate

# Directory where the examples in Literate.jl format are
lectures_DIR = abspath(@__DIR__, "lectures")
lectures = [f for f in readdir(lectures_DIR) if startswith(f, "lecture")]
# Directory where the examples will be generated
generated_DIR = abspath(joinpath(@__DIR__, "generated"))
# Erase previous generated contents if they exist
isdir(generated_DIR) && rm(generated_DIR, recursive=true, force=true)
# Create generated if it did not exist
isdir(generated_DIR) || mkdir(generated_DIR)

println("Generating files using Literate.jl")
for lecture in lectures
    println("-----------------------------")
    println("$lecture")
    for f in readdir(joinpath(lectures_DIR, lecture))
        input = joinpath(lectures_DIR, lecture, f)
        output_DIR = joinpath(generated_DIR, lecture)
        mkdir(output_DIR) # create lecture DIR
        if endswith(f, ".jl")
            println("\ncreating $f:")
            # Alternate the commented line to avoid executing the notebook
            # when you are editing the text only
            #Literate.notebook(input, output_DIR, execute = true)
            Literate.notebook(input, output_DIR, execute = false)
        elseif endswith(f, ".md")
            println("\ncopying $f")
            cp(input, joinpath(output_DIR, f))
        end
    end
end

