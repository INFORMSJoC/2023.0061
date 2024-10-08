using BSTPDDD
using Test

@testset "BSTPDDD.jl" begin
    # test solution for instance I01
    rsl = BSTPDDD.ddd_bstp("I01",15,5);
    @test rsl[4] == 270;
end