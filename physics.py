import math
import numpy as np


class Physics:
    def calculate_buoyancy(self, vol, densflu):
        """
        Calculates the buoyancy force experienced by an object submerged in a fluid.

        Args:
            vol (float): Volume of the object in cubic meters.
            densflu (float): Density of the fluid in kilograms per cubic meter.

        Returns:
            float: Buoyancy force in newtons.
        """
        g = 9.81  # m/s ** 2
        buoyancy = g * vol * densflu
        return buoyancy

    def will_it_float(self, vol, mass):
        """
        Determines whether an object will float or sink in water based on its volume and mass.

        Args:
            vol (float): Volume of the object in cubic meters.
            mass (float): Mass of the object in kilograms.

        Returns:
            bool: True if the object will float, False if it will sink.
        """
        g = 9.81  # m/s ** 2
        density_water = 1000  # kg/m ** 3 (density of water)

        buoyancy_force = self.calculate_buoyancy(vol, density_water)
        gravitational_force = mass * g

        if buoyancy_force >= gravitational_force:
            return True
        else:
            return False

    def calculate_pressure(self, depth):
        """
        Calculates the pressure exerted by a column of water at a given depth.

        Args:
            depth (float): Depth of the water column in meters.

        Returns:
            float: Pressure exerted by the water column in pascals.
        """
        g = 9.81  # m/s ** 2
        density_water = 1000  # kg/m ** 3 (density of water)
        pressure = density_water * g * depth
        return pressure

    def calculate_acceleration(self, force, mass):
        """
        Calculates the acceleration of an object given the force applied to it and its mass.

        Args:
            F (float): Force applied to the object in Newtons.
            m (float): Mass of the object in kilograms.

        Returns:
            float: Acceleration of the object in meters per second squared (m/s^2).
        """
        acceleration = force / mass
        return acceleration

    def calculate_angular_acceleration(self, inertia, torque):
        """
        Calculates the angular acceleration of an object given the torque applied to it and its moment of inertia.

        Args:
            tau (float): Torque applied to the object in Newton-meters.
            I (float): Moment of inertia of the object in kilogram-meters squared (kg⋅m^2).

        Returns:
            float: Angular acceleration of the object in radians per second squared (rad/s^2).
        """
        torque = inertia * angular_acceleration
        angular_acceleration = torque / inertia
        return angular_acceleration

    def calculate_torque(self, r, force):
        """
        Calculates the torque exerted on an object given the distance from the pivot point and the applied force.

        Args:
            r (float): Distance from the pivot point (lever arm) in meters.
            force (float): Applied force in Newtons.

        Returns:
            float: Torque exerted on the object in Newton-meters.
        """
        torque = r * force
        return torque

    def calculate_moment_of_inertia(self, mass, r):
        """
        Calculates the moment of inertia of an object given its mass and the distance from the axis of rotation.

        Args:
            mass (float): Mass of the object in kilograms.
            r (float): Distance from the axis of rotation in meters.

        Returns:
            float: Moment of inertia of the object in kilogram-meters squared (kg⋅m^2).
        """
        moment_of_inertia = mass * r**2
        return moment_of_inertia

    def calculate_auv_acceleration(
        self, F_magnitude, F_angle, mass=100, volume=0.1, thruster_distance=0.5
    ):
        """
        Calculates the acceleration of the AUV in the 2D plane.

        Args:
            F_magnitude (float): Magnitude of force applied by the thruster in Newtons.
            F_angle (float): Angle of the force applied by the thruster in radians.
                            Positive angles are measured in the counter-clockwise direction from the x-axis.
            mass (float, optional): Mass of the AUV in kilograms. Defaults to 100 kg.
            volume (float, optional): Volume of the AUV in cubic meters. Defaults to 0.1 m^3.
            thruster_distance (float, optional): Distance from the center of mass of the AUV to the thruster in meters.
                                                Defaults to 0.5 m.

        Returns:
            float: Acceleration of the AUV in meters per second squared.
        """
        density_water = 1000  # kg/m^3 (density of water)
        g = 9.81  # m/s^2 (acceleration due to gravity)

        # Calculate buoyancy force
        buoyancy_force = density_water * g * volume

        # Calculate gravitational force
        gravitational_force = mass * g

        # Calculate net force acting on the AUV
        net_force = F_magnitude * math.cos(F_angle)

        # Calculate acceleration using Newton's second law (F = m*a)
        acceleration = net_force / mass

        return acceleration

    def calculate_auv2_acceleration(self, T, alpha, theta, mass=100):
        """
        Calculates the acceleration of the AUV in the 2D plane.

        Args:
            T (np.ndarray): Magnitudes of the forces applied by the thrusters in Newtons.
            alpha (np.ndarray): Angles of the thrusters in radians.
            theta (float): Angle of the AUV in radians.
            mass (float, optional): Mass of the AUV in kilograms. Default value is 100 kg.

        Returns:
            np.ndarray: Acceleration of the AUV in the 2D plane as a numpy array [ax, ay].
        """
        # Calculate the total force components in the AUV's frame
        Fx_prime = np.array(
            [np.cos(alpha), np.cos(alpha), -1 * np.cos(alpha), -1 * np.cos(alpha)]
        )
        Fy_prime = np.array(
            [np.sin(alpha), -1 + np.sin(alpha), -1 * np.sin(alpha), np.sin(alpha)]
        )

        # Calculate the rotational matrix from the AUV's frame to the world frame
        R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

        # Calculate the total force components in the world frame
        F_vec = np.array([Fx_prime, Fy_prime])
        F_world = np.dot(R, F_vec)

        # Calculate the acceleration components using Newton's second law (F = m*a)
        acceleration = F_world / mass

        return acceleration

    def calculate_auv2_angular_acceleration(self, T, alpha, L, inertia=100):
        """
        Calculates the angular acceleration of the AUV.

        Args:
            T (np.ndarray): Magnitudes of the forces applied by the thrusters in Newtons.
            alpha (np.ndarray): Angles of the thrusters in radians.
            L (float): Distance from the center of mass of the AUV to the thrusters in meters.
            l (float): Distance from the center of mass of the AUV to the thrusters in meters.
            inertia (float, optional): Moment of inertia of the AUV in kg * m^2. Default value is 100.

        Returns:
            float: Angular acceleration of the AUV in radians per second squared (rad/s^2).
        """
        # Calculate the torques produced by each thruster
        torque = T * L * np.sin(alpha)

        # Calculate the net torque
        net_torque = np.sum(torque)

        # Calculate the angular acceleration
        angular_acceleration = net_torque / inertia

        return angular_acceleration


var = Physics()
print(var.calculate_buoyancy(30, 25))
print(var.will_it_float(25, 15))
print(var.calculate_pressure(50))
print(var.calculate_acceleration(45, 35))
print(var.calculate_angular_acceleration(15, 20))
print(var.calculate_torque(15, 50))
print(var.calculate_moment_of_inertia(250, 10))
print(var.calculate_auv_acceleration(100, 45, 100, 0.1, 0.5))
print(var.calculate_auv2_acceleration(55, 35, 45, 100))
print(var.calculate_auv2_angular_acceleration(25, 35, 15, 100))
