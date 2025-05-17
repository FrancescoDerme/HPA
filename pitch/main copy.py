from manim import *

class scene1(Scene):
    def construct(self):

        title_text = Tex("hp-adaptive strategies for high-fidelity multi-physics simulations", font_size=35, color=WHITE)
        self.play(Write(title_text), run_time=3) # Corresponds to "we'll present our plan to build a"
        self.wait(0.8) # Hold the title for a moment

        # --- Author Names (0.5s - 3.0s) ---
        authors_name = Text("Pietro Fumagalli       Francesco Derme", font_size=25, color=BLUE_C)
        authors_group = VGroup(authors_name).arrange(RIGHT, buff=0.25)
        authors_group.move_to(DOWN*0.7) # Positioned below the title

        self.play(Write(authors_name), run_time=2) # Corresponds to "I'm Pietro Fumagalli"
        self.wait(0.2) # slight pause before Francesco's name

        # --- Transition to Plan (3.0s - 3.5s) ---
        # Authors and hello shift up to make space
        main_titles_group = VGroup(title_text, authors_group)

        self.play(FadeOut(main_titles_group, shift=UP*0.5), run_time=1.0)
        self.wait(0.2)

        # 2. Domain Boundary Appears (Using a hexagon as the domain Omega)
        domain_boundary = RegularPolygon(n=6, radius=1.7, color=BLUE_D, stroke_width=4)
        domain_boundary.move_to(ORIGIN)
        self.play(Create(domain_boundary), run_time=1.0)
        self.wait(0.5)

        # 3. Coarse Mesh Lines Appear
        # These lines will divide the hexagon into 6 large triangular elements.
        center = domain_boundary.get_center()
        vertices = domain_boundary.get_vertices()
        
        coarse_mesh_lines = VGroup()
        for vertex in vertices:
            coarse_mesh_lines.add(Line(center, vertex, color=WHITE, stroke_width=2.5))
        
        self.play(Create(coarse_mesh_lines), run_time=1.0)
        self.wait(0.5)

        # 4. Fine Mesh Lines Appear (Refinement of the coarse triangular elements)
        # Each large triangle will be subdivided into 4 smaller triangles.
        fine_mesh_lines = VGroup()
        for i in range(6):
            # Vertices of the current coarse triangle: center, vertices[i], vertices[(i+1)%6]
            v0 = center
            v1 = vertices[i]
            v2 = vertices[(i + 1) % 6] # Ensures the last vertex connects to the first

            # Midpoints of the edges of this coarse triangle
            m01 = (v0 + v1) / 2  # Midpoint of line from center to vertex i
            m12 = (v1 + v2) / 2  # Midpoint of outer edge
            m02 = (v0 + v2) / 2  # Midpoint of line from center to vertex i+1

            # Lines connecting these midpoints form the refinement
            fine_mesh_lines.add(Line(m01, m12, color=YELLOW_C, stroke_width=1.5))
            fine_mesh_lines.add(Line(m12, m02, color=YELLOW_C, stroke_width=1.5))
            fine_mesh_lines.add(Line(m02, m01, color=YELLOW_C, stroke_width=1.5))
        
        # Animate the appearance of these finer mesh lines
        self.play(LaggedStart(*[Create(line) for line in fine_mesh_lines], lag_ratio=0.05), run_time=2.0)
        self.wait(0.5)
        
        # Group all mesh lines for the fade-out animation
        all_mesh_elements = VGroup(domain_boundary, coarse_mesh_lines, fine_mesh_lines)
        
        self.play(FadeOut(all_mesh_elements, scale=3.5), run_time=1.5)


import numpy as np

class scene3(ThreeDScene):
    def construct(self):
        # === Part 1: Title, Equation, and Boundary Description (2D Camera View) ===
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES)
        self.camera.frame_height = 8.5

        # 1. Title
        title = Tex("A first example: the Laplacian", font_size=46)
        title.move_to(ORIGIN)
        self.play(FadeIn(title, scale=0.7), run_time=1.0)
        self.play(title.animate.scale(0.9).to_edge(UP, buff=0.4), run_time=1.0)
        self.wait(0.2)

        # 2. Laplace Equation System (No Neumann condition)
        laplace_eq_system = MathTex(
            r"""\begin{cases}
                -\Delta u(\mathbf{x}) = f(\mathbf{x}) & \mathbf{x} \in \Omega \\
                u(\mathbf{x}) = g(\mathbf{x}) & \mathbf{x} \in \partial\Omega 
            \end{cases}""", # Using dOmega for general Dirichlet BC
            font_size=38
        ).next_to(title, DOWN, buff=0.35)
        self.play(Write(laplace_eq_system), run_time=2.5) # Adjusted runtime
        self.wait(0.2)

        # 3. Explicit Boundary Expression for a Square Domain (Centered)
        L_char = "L"
        boundary_expr = MathTex(
            r"""\text{For } \Omega = (-""" + L_char + r""", """ + L_char + r""")^2 \text{, the boundary } \partial\Omega = \bigcup_{i=1}^4 \Gamma_i \text{ is:}""",
            font_size=30 # Smaller font for this detailed part
        )
        # Position it centered below the Laplace equation system
        boundary_expr.next_to(laplace_eq_system, DOWN, buff=0.5)

        boundary_details = MathTex(            
            r"""\begin{alignedat}{2}
                \Gamma_1 &= \{ (""" + L_char + r""", y)   &&\mid -""" + L_char + r""" \le y \le """ + L_char + r""" \} \\
                \Gamma_2 &= \{ (-""" + L_char + r""", y)  &&\mid -""" + L_char + r""" \le y \le """ + L_char + r""" \} \\
                \Gamma_3 &= \{ (x, """ + L_char + r""")   &&\mid -""" + L_char + r""" \le x \le """ + L_char + r""" \} \\
                \Gamma_4 &= \{ (x, -""" + L_char + r""")  &&\mid -""" + L_char + r""" \le x \le """ + L_char + r""" \}
            \end{alignedat}""",
            font_size=30
        )
     
        self.play(Write(boundary_expr), run_time=4.0)
        boundary_details.next_to(boundary_expr, DOWN, buff=0.5)
        self.play(Write(boundary_details), run_time=4.0)
        self.wait(1.0)

        # 4. Disappearance of all text elements
        text_elements = VGroup(title, laplace_eq_system, boundary_expr, boundary_details)
        self.play(FadeOut(text_elements, shift=UP * 0.3), run_time=1.0)

        # === Part 2: 3D Plot with Non-Homogeneous Mesh & Modified Solution ===
        self.move_camera(phi=60 * DEGREES, theta=-55 * DEGREES, distance=14, run_time=1.8)

        L_val = 1.0 # Numerical value for L

        # 5a. 3D Cartesian Grid
        axes = ThreeDAxes(
            x_range=[-L_val*1.2, L_val*1.2, L_val/2], 
            y_range=[-L_val*1.2, L_val*1.2, L_val/2], 
            z_range=[-0.3, 1.3, 0.2], # Adjusted for potential noise amplitude
            x_length=5.5, y_length=5.5, z_length=3.8,
            axis_config={"include_numbers": False, "font_size": 18, "include_tip": False},
        )
        axes_labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u")
        self.play(Create(axes), Write(axes_labels), run_time=1.5)

        # 5b. Boundary of the square on XY-plane
        domain_boundary_3d = Polygon(
            axes.c2p(-L_val, -L_val, 0), axes.c2p(L_val, -L_val, 0),
            axes.c2p(L_val, L_val, 0), axes.c2p(-L_val, L_val, 0),
            color=YELLOW_C, stroke_width=4, fill_opacity=0
        )
        self.play(Create(domain_boundary_3d), run_time=1.0)
        self.wait(0.3)

        # 5c. Non-homogeneous mesh lines
        mesh_lines = VGroup()
        mesh_lines.add(Line(axes.c2p(-L_val, 0, 0), axes.c2p(L_val, 0, 0), stroke_color=GREY_B, stroke_width=1.5)) # x-axis line
        mesh_lines.add(Line(axes.c2p(0, -L_val, 0), axes.c2p(0, L_val, 0), stroke_color=GREY_B, stroke_width=1.5)) # y-axis line

        fine_step = L_val / 8.0
        coarse_step = L_val / 2.0

        # Bottom-left quadrant (x from -L to 0, y from -L to 0) -> FINE MESH
        for i in range(1, int(L_val / fine_step)): # From 1 to exclude boundary lines already in domain_boundary_3d or central lines
            x_coord_bl = -L_val + i * fine_step
            y_coord_bl = -L_val + i * fine_step
            if x_coord_bl < 0: # strictly within the -L to 0 part for x
                 mesh_lines.add(Line(axes.c2p(x_coord_bl, -L_val, 0), axes.c2p(x_coord_bl, 0, 0), stroke_color=WHITE, stroke_width=1.2))
            if y_coord_bl < 0: # strictly within the -L to 0 part for y
                 mesh_lines.add(Line(axes.c2p(-L_val, y_coord_bl, 0), axes.c2p(0, y_coord_bl, 0), stroke_color=WHITE, stroke_width=1.2))
        
        # Other quadrants (coarser) - example for Top-Right (x from 0 to L, y from 0 to L)
        for i in range(1, int(L_val / coarse_step)):
            x_coord_tr = 0 + i * coarse_step
            y_coord_tr = 0 + i * coarse_step
            if x_coord_tr < L_val:
                mesh_lines.add(Line(axes.c2p(x_coord_tr, 0, 0), axes.c2p(x_coord_tr, L_val, 0), stroke_color=WHITE, stroke_width=1.2))
            if y_coord_tr < L_val:
                mesh_lines.add(Line(axes.c2p(0, y_coord_tr, 0), axes.c2p(L_val, y_coord_tr, 0), stroke_color=WHITE, stroke_width=1.2))
        # Note: You can fill in the other two quadrants (bottom-right, top-left) with coarse lines similarly.
        # This example just shows one fine and one coarse quadrant explicitly defined beyond the central lines.

        self.play(Create(mesh_lines, lag_ratio=0.1), run_time=2.5)
        self.wait(0.5)

        # 5d. Solution Plot with modified behavior
        def piecewise_solution_func(x, y):
            peak_amplitude = 0.9 # Max height of peaks
            sigma_peak = L_val * 0.20 # Sharpness of peaks
            oscillation_freq = 18 / L_val # For "problematic" oscillations
            oscillation_amplitude_factor = 0.3 # Factor for oscillation height relative to peak
            constant_value_coarse = 0.1 # Value for constant regions

            val = constant_value_coarse # Default for coarse regions (TL, BR)

            # Bottom-Left: "Problematic" (peak with oscillations)
            if x <= 0.001 and y <= 0.001: # (x <= 0 and y <= 0) with tolerance
                center_x, center_y = -L_val * 0.6, -L_val * 0.6
                dist_sq = (x - center_x)**2 + (y - center_y)**2
                gaussian_shape = np.exp(-dist_sq / (2 * sigma_peak**2))
                oscillations = np.sin(oscillation_freq * x) * np.cos(oscillation_freq * y)
                val = peak_amplitude * gaussian_shape * (1 + oscillation_amplitude_factor * oscillations)
            
            # Top-Right: "Problematic" (peak with different oscillations)
            elif x > -0.001 and y > -0.001: # (x > 0 and y > 0) with tolerance
                center_x, center_y = L_val * 0.6, L_val * 0.6
                dist_sq = (x - center_x)**2 + (y - center_y)**2
                gaussian_shape = np.exp(-dist_sq / (2 * sigma_peak**2))
                oscillations = np.cos(oscillation_freq * x) * np.sin(oscillation_freq * y * 1.2) # Slightly different pattern
                val = peak_amplitude * gaussian_shape * (1 + oscillation_amplitude_factor * oscillations)
            
            # TL (x<=0, y>0) and BR (x>0, y<=0) regions will use constant_value_coarse
            
            return val

        solution_surface = Surface(
            lambda u, v: axes.c2p(u, v, piecewise_solution_func(u, v)),
            u_range=[-L_val, L_val], v_range=[-L_val, L_val],
            resolution=(72, 72), # Higher resolution for detail
            fill_opacity=0.85,
            # Colors that might show detail well
            # Using a colormap can also be an option via set_shader_input
            checkerboard_colors=[BLUE_E, BLUE_C], 
            stroke_color=BLACK,
            stroke_width=0.1
        )
        self.play(Create(solution_surface), run_time=4.0) # Longer time for complex surface
        self.wait(4.0) # Hold the final plot longer

        # No specific text label for problematic regions as per request.
        # The visual difference in the surface itself is the cue.

        # Cleanup
        all_3d_elements = VGroup(axes, axes_labels, domain_boundary_3d, mesh_lines, solution_surface)
        self.play(FadeOut(all_3d_elements), run_time=1.5)
        self.wait(0.5)