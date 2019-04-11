# MMPDElab
MMPDElab package
MMPDElab is a package written in MATLAB for adaptive mesh movement and adaptive moving mesh
    P1 finite element solution of second-order partial different equations (PDEs) having continuous
    solutions. The adaptive mesh movement is based on the new implementation of the moving mesh partial
    differential equation (MMPDE) method. The mesh equation is integrated using either ode45
    (an explicit MATLAB ODE solver) or ode15s (an implicit MATLAB ODE solver) while physical PDEs
    are discretized in space using P1 conforming finite elements on moving meshes and integrated
    in time with the fifth-order Radau IIA method (an implicit Runge-Kutta method) with a two-step
    error estimator for time step selection. An introduction of this package is given
    in intro_MMPDElab.pdf (arXiv:XXXXX) contained in this distribution.
    
    MMPDElab is a package written in MATLAB  for adaptive mesh 
    movement and adaptive moving mesh P1 finite element solution 
    of partial different equations having continuous solutions.
    
    Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    MMPDElab is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    MMPDElab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU Affero General Public License at
    <https://www.gnu.org/licenses/>.
