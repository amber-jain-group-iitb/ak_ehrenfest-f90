program ehrenfest_dynamics

    ! generate pes
    use model_spin_boson
    call pes(1)

    ! !! run dynamics
    ! use dynamics    ! specify model to use in mod_dynamics_ef 
    ! call setup()
    ! call run_dynamics()
    
end program ehrenfest_dynamics
