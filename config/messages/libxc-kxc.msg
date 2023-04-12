The LibXC package does not provide the 3rd derivatives of
energy!
You should recompile it using --enable-kxc configure option.
Otherwise Abinit will not be able to provide nonlinear
response properties \\(Raman intensities, ...\\).
