using HubbardRenyiDQMC

# Sites are numbered x + (y - 1)*Lx. Change the size and boundary conditions here.
K = square_lattice_kinetic_matrix(2, 2; t=1.0, periodic=(false, false))
ground = ed_ground_state(K; U=4.0)
thermal = ed_thermal(K; U=4.0, beta=2.0)

println("Ground-state energy per site: ", ground.energy_per_site)
println("Left-column Renyi entropy S2: ", ed_renyi2(ground, [1, 3]))
println("Thermal energy per site at beta=2: ", thermal.energy_per_site)
println("Thermal density: ", thermal.density)
println("Thermal double occupancy per site: ", thermal.double_occupancy)
