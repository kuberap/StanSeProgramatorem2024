import numpy as np
import vtk
import matplotlib.pyplot as plt

import PIL.Image as Image
class LBM:
    # Private fix parameters
    __p0 = 0
    __nu_LB = 0.01
    __cs2 = 1. / 3.
    __tau_LB = (__nu_LB + 0.5) / __cs2



    def __init__(self, uMax, nu_phys, Lx, Ly, T, dx, Cx = 0.2, Cy = 0.2, R=0.05):
        """
        Initialize the LBM simulation.
        """
        self.Lx = Lx
        self.Ly = Ly
        self.uMax = uMax
        self.dx = dx

        # cylinder
        self.Cx = Cx
        self.Cy = Cy
        self.R = R
        #-------------------------------------------

        self.nx = int(np.ceil(Lx / dx))
        self.ny = int(np.ceil(Ly / dx))

        # diffusive scaling
        self.dt = self.__nu_LB / nu_phys * dx * dx
        self.currentIteration = 0
        self.numberOfItereations = int(np.ceil(T / self.dt))

        # tady jsem to musel priohnout
        self.nu = nu_phys  # nunto pridat
        self.u0 = self.phys2LBMVelocity(self.uMax)
        # -------
        print(f"uMAx: {self.uMax}, nu: {self.nu}, Lx: {Lx}, Ly: {Ly}, T: {T}, dx: {dx}")
        print(f"nx: {self.nx}, ny: {self.ny}, dt: {self.dt}, numberOfItereations: {self.numberOfItereations}")
        #exit(0)
        # macroscopic variables
        self.density = np.array(np.zeros((self.nx * self.ny), dtype=float))
        self.velocity = np.array(np.zeros((self.nx * self.ny, 2), dtype=float))
        self.WallMap = np.array(np.zeros((self.nx * self.ny), dtype=int))

        # main variables
        self.df1 = np.array(np.zeros((self.nx * self.ny, 9), dtype=float))
        self.df2 = np.array(np.zeros((self.nx * self.ny, 9), dtype=float))

        # mapping for wallMap
        self.flowID = 0
        self.wallID = 1
        self.inflowID = 2
        self.outflowID = 3

        # Model parameters
        # relaxation frequency
        self.omega = (self.__cs2) / (self.__nu_LB + 0.5)  # Relaxation parameter
        self.weights = [4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.]
        self.velocitySet = np.transpose(
            [[0., 0.], [1., 0.], [0., 1.], [-1., 0.], [0., -1.], [1., 1.], [-1., 1.], [-1., -1.], [1., -1.]])

    # Functions for collision
    def equilibrium(self, i_rho, i_u, i_v):
        eq = []
        for idx in range(9):
            eq.append(
                i_rho * self.weights[idx] * (1. + 3 * (self.velocitySet[0][idx] * i_u + self.velocitySet[1][idx] * i_v)
                                             + 9. / 2. * (self.velocitySet[0][idx] * i_u + self.velocitySet[1][
                            idx] * i_v) ** 2
                                             - 3. / 2. * (i_u ** 2 + i_v ** 2)))
        return np.array(eq)

    # compute macroscopic variables
    def computeMacro(self, i_df):
        rho_loc = np.sum(i_df)
        vx_loc = np.sum(i_df * self.velocitySet[0]) / rho_loc
        vy_loc = np.sum(i_df * self.velocitySet[1]) / rho_loc
        return rho_loc, np.array([vx_loc, vy_loc])

    # apply boundary conditions
    # inflow boundary condition
    def inflowBoundary(self, i_x, i_y, i_df):
        xf, yf = self.positionCellCenter(i_x, i_y)
        loc_df_plus = self.streaming(i_df, i_x + 1, i_y)
        rho_plus, _ = self.computeMacro(loc_df_plus)
        no1Ocs2 = 3.
        InputFactor = no1Ocs2 * self.phys2LBMVelocity(self.phys2LBMVelocity(self.nu * 8. * self.u0 / self.Ly ** 2))
        phys_vx, phys_vy = self.inflowVelocity(yf)
        loc_vx = self.phys2LBMVelocity(phys_vx)
        loc_vy = self.phys2LBMVelocity(phys_vy)
        loc_rho = rho_plus - InputFactor
        return self.equilibrium(loc_rho, loc_vx, loc_vy), loc_rho, np.array([loc_vx, loc_vy])

    # bounce back boundary condition
    def bounceback(self, i_df):
        out_df = [i_df[0]]
        out_df.append(i_df[3])
        out_df.append(i_df[4])
        out_df.append(i_df[1])
        out_df.append(i_df[2])
        out_df.append(i_df[7])
        out_df.append(i_df[8])
        out_df.append(i_df[5])
        out_df.append(i_df[6])
        return out_df

    # outflow boundary condition
    def outflowCondition(self, i_x, i_y, i_df):
        return i_df

    # Streaming
    def streaming(self, i_df, i_x, i_y):
        out_df = [i_df[self.D1Idx(i_x, i_y)][0]]
        xm = np.max([0, (i_x - 1)])
        ym = np.max([0, (i_y - 1)])
        xp = np.min([self.nx - 1, (i_x + 1)])
        yp = np.min([self.ny - 1, (i_y + 1)])
        out_df.append(i_df[self.D1Idx(xm, i_y)][1])
        out_df.append(i_df[self.D1Idx(i_x, ym)][2])
        out_df.append(i_df[self.D1Idx(xp, i_y)][3])
        out_df.append(i_df[self.D1Idx(i_x, yp)][4])
        out_df.append(i_df[self.D1Idx(xm, ym)][5])
        out_df.append(i_df[self.D1Idx(xp, ym)][6])
        out_df.append(i_df[self.D1Idx(xp, yp)][7])
        out_df.append(i_df[self.D1Idx(xm, yp)][8])
        return out_df

    # Update step
    def updateStep(self, i_df, i_x, i_y):
        # streaming
        loc_df = self.streaming(i_df, i_x, i_y)
        mapID = self.WallMap[self.D1Idx(i_x, i_y)]
        # perform boundary condition
        if mapID == self.wallID:
            loc_df = self.bounceback(loc_df)
            loc_density = 1
            loc_velocity = np.array([0., 0.])
        # perform initial condition
        elif mapID == self.inflowID:
            loc_df, loc_density, loc_velocity = self.inflowBoundary(i_x, i_y, i_df)
        # perform outflow condition
        elif mapID == self.outflowID:
            loc_df = self.outflowCondition(i_x, i_y, loc_df)
            loc_density, loc_velocity = self.computeMacro(loc_df)
        else:
            # compute macro
            loc_density, loc_velocity = self.computeMacro(loc_df)
            # compute collision
            loc_df = loc_df + self.omega * (self.equilibrium(loc_density, loc_velocity[0], loc_velocity[1]) - loc_df)
            # return i_df, i_density, i_velocity
        return loc_df, loc_density, loc_velocity

    # Go to next iteration
    def nextIteration(self):
        t = self.currentIteration
        if (t < self.numberOfItereations):
            self.currentIteration = self.currentIteration + 1
            for x in range(self.nx):
                for y in range(self.ny):
                    idx = self.D1Idx(x, y)
                    if (t % 2 == 0):
                        self.df2[idx], self.density[idx], self.velocity[idx] = self.updateStep(self.df1, x, y)
                    else:
                        self.df1[idx], self.density[idx], self.velocity[idx] = self.updateStep(self.df2, x, y)

    # get physical time
    def physicalTime(self):
        return self.currentIteration * self.dt

    # mapping to 1D array
    def D1Idx(self, i_x, i_y):
        return i_x + i_y * self.nx

    # return position of cell center
    def positionCellCenter(self, i_x, i_y):
        return i_x * self.dx, i_y * self.dx - self.dx / 2.

    # return position of cell bottom left corner
    def positionCell(self, i_x, i_y):
        return i_x * self.dx - self.dx / 2., i_y * self.dx - self.dx

    # scaling to non-dimensional units
    def phys2LBMVelocity(self, i_u):
        return i_u * self.dt / self.dx

    # scaling to physical units
    def LBM2PhysVelocity(self, i_u):
        return i_u * self.dx / self.dt

    # initialize map
    def mapInit(self):
        for x in range(self.nx):
            idxTop = self.D1Idx(x, self.ny - 1)
            idxDown = self.D1Idx(x, 0)
            self.WallMap[idxTop] = self.wallID
            self.WallMap[idxDown] = self.wallID

        for y in range(1, self.ny - 1):
            idxLeft = self.D1Idx(0, y)
            idxRight = self.D1Idx(self.nx - 1, y)
            self.WallMap[idxLeft] = self.inflowID
            self.WallMap[idxRight] = self.outflowID

        for x in range(self.nx):
            for y in range(self.ny):
                xf, yf = self.positionCellCenter(x, y)
                if ((xf - self.Cx) ** 2 + (yf - self.Cy) ** 2 <= self.R ** 2): # zmena
                    idx1d = self.D1Idx(x, y)
                    self.WallMap[idx1d] = self.wallID

    # initialize macroscopic variables
    def macroInit(self):
        # initialization
        for x in range(self.nx):
            for y in range(self.ny):
                idx1d = self.D1Idx(x, y)
                self.velocity[idx1d][0] = self.phys2LBMVelocity(0.)
                self.velocity[idx1d][1] = self.phys2LBMVelocity(0.)
                self.density[idx1d] = 1.

    # initialize distribution functions
    def distributionInit(self):
        # initial condition
        for x in range(self.nx):
            for y in range(self.ny):
                idx1d = self.D1Idx(x, y)
                self.df1[idx1d] = self.equilibrium(self.density[idx1d], self.velocity[idx1d][0],
                                                   self.velocity[idx1d][1])

    # start simulation
    def start(self):
        self.mapInit()
        self.macroInit()
        self.distributionInit()

    # inflow velocity profile
    def inflowVelocity(self, i_y):
        u0 = self.phys2LBMVelocity(self.uMax)
        return 4. * u0 * i_y * (self.Ly - i_y) / self.Ly ** 2, 0

    # compute macroscopic variables from distribution function in physical units
    def computePhysMacro(self, i_df):
        nonDimCS2 = self.phys2LBMVelocity(self.phys2LBMVelocity(1. / 3.))
        rho_loc = np.sum(i_df)
        vx_loc = np.sum(i_df * self.velocitySet[0]) / rho_loc
        vy_loc = np.sum(i_df * self.velocitySet[1]) / rho_loc
        return (rho_loc - 1) * self.__cs2 + self.__p0, np.array(
            [self.LBM2PhysVelocity(vx_loc), self.LBM2PhysVelocity(vy_loc)])

        # write vtk output

    def writeVTK(self, i_name):
        # Create the rectilinear grid
        rect_grid = vtk.vtkRectilinearGrid()

        # Define the number of points along each axis
        x_coords = vtk.vtkFloatArray()
        y_coords = vtk.vtkFloatArray()
        z_coords = vtk.vtkFloatArray()

        # Define the coordinates
        x_coords.SetNumberOfValues(self.nx + 1)
        y_coords.SetNumberOfValues(self.ny + 1)
        z_coords.SetNumberOfValues(1)

        for x in range(self.nx + 1):
            x_coords.SetValue(x, self.positionCell(x, 0)[0])
        for y in range(self.ny + 1):
            y_coords.SetValue(y, self.positionCell(0, y)[1])

        z_coords.SetValue(0, 0.0)

        # Set the coordinates in the rectilinear grid
        rect_grid.SetDimensions((self.nx + 1), (self.ny + 1), 1)
        rect_grid.SetXCoordinates(x_coords)
        rect_grid.SetYCoordinates(y_coords)
        rect_grid.SetZCoordinates(z_coords)

        # Create the scalar array
        mapArray = vtk.vtkIntArray()
        mapArray.SetNumberOfValues(self.nx * self.ny)  # 4 * 3 * 2 = 24

        # Create the scalar array
        pressureArray = vtk.vtkFloatArray()
        pressureArray.SetNumberOfValues(self.nx * self.ny)  # 4 * 3 * 2 = 24

        # pressureArray_an = vtk.vtkFloatArray()
        # pressureArray_an.SetNumberOfValues(i_Nx*i_Ny)  # 4 * 3 * 2 = 24

        velocity = vtk.vtkFloatArray()
        velocity.SetNumberOfComponents(3)  # 3 components for 3D vectors (x, y, z)
        velocity.SetNumberOfTuples(self.nx * self.ny)  # Same number of points as the grid

        # velocity_an = vtk.vtkFloatArray()
        # velocity_an.SetNumberOfComponents(3)  # 3 components for 3D vectors (x, y, z)
        # velocity_an.SetNumberOfTuples(i_Nx*i_Ny)  # Same number of points as the grid

        # Fill the scalar array with some values
        nonDimCS2 = self.phys2LBMVelocity(self.phys2LBMVelocity(1. / 3.))
        for x in range(self.nx):
            for y in range(self.ny):
                mapArray.SetValue(x + y * self.nx, self.WallMap[self.D1Idx(x, y)])
                pressureArray.SetValue(x + y * self.nx, (self.density[self.D1Idx(x, y)] - 1.) * nonDimCS2)
                velocity.SetTuple3(x + y * self.nx, self.LBM2PhysVelocity(self.velocity[self.D1Idx(x, y)][0]),
                                   self.LBM2PhysVelocity(self.velocity[self.D1Idx(x, y)][1]), 0)

        # Attach the scalar array to the rectilinear grid
        mapArray.SetName("Map")
        rect_grid.GetCellData().AddArray(mapArray)

        # Attach the scalar array to the rectilinear grid
        pressureArray.SetName("Pressure")
        rect_grid.GetCellData().AddArray(pressureArray)

        # Attach the vector array to the rectilinear grid
        velocity.SetName("velocity")
        rect_grid.GetCellData().AddArray(velocity)

        # Write the rectilinear grid to a VTK file
        writer = vtk.vtkRectilinearGridWriter()
        writer.SetFileName(i_name)
        writer.SetInputData(rect_grid)
        writer.Write()

    # nunto udelat property pro výstup
    @property
    def velocity2D(self):
        return self.velocity.reshape((self.ny,self.nx, 2)).copy()

    @property
    def velocity2DImage(self):
        return Image.fromarray(np.sqrt(self.velocity2D[:,:,0]**2+self.velocity2D[:,:,1]**2))

if __name__ == '__main__':
    # Parameters
    Lx = 2.2
    Ly = 0.41
    T = 1000  # maximal time
    nu = 0.001  # for test case only, no physical meening
    u0 = 0.3  # 1m/s
    rho_0 = 1
    dx = 0.01
    # cylinder
    Cx = 0.2
    Cy = 0.2
    R = 0.05

    # Create LBM object
    instance = LBM(u0,nu, Lx, Ly,T, dx)
    instance.start()
    print(instance.velocity2D.shape)
    while instance.physicalTime()<T:
        instance.nextIteration()
        print(instance.physicalTime())
        plt.imshow(np.sqrt(instance.velocity2D[:,:,0]**2+instance.velocity2D[:,:,1]**2))
        plt.show()

