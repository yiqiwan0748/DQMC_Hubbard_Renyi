using Test

# Evaluate the actual reservation helper without loading optional physics libraries.
source = joinpath(@__DIR__, "..", "..", "examples", "dmrg_square.jl")
expressions = Meta.parseall(read(source, String)).args
helper = only(filter(expressions) do expression
    expression isa Expr && expression.head == :function &&
        expression.args[1] isa Expr &&
        expression.args[1].args[1] == :reserve_dmrg_output
end)
eval(helper)

println("Julia ", VERSION, "; kernel ", Sys.KERNEL)
@testset "Native output-file reservation" begin
    mktempdir() do directory
        output = joinpath(directory, "nested", "result.csv")
        io = reserve_dmrg_output(output)
        @test isfile(output)
        @test_throws ErrorException reserve_dmrg_output(output)
        write(io, "sentinel\n")
        close(io)
        @test read(output, String) == "sentinel\n"
        @test_throws ErrorException reserve_dmrg_output(output)
        @test read(output, String) == "sentinel\n"
    end
end
