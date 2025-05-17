from manim import *
import numpy as np
import random
import math

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

class scena2(Scene):
    def construct(self):
        h = Text("h", color=BLUE_C, font_size = 45).to_corner(UP)
        p = Text("p", color=RED_C, font_size = 45).to_corner(UP)
        x_position = config.frame_width / 4
        h.move_to([-x_position, h.get_y(), 0])
        p.move_to([x_position, p.get_y(), 0])

        start_pos = [0, h.get_y(), 0]
        end_pos = [0, -h.get_y(), 0]
        middle_line = Line(start=start_pos, end=end_pos)

        self.play(Write(h),
                  Write(p),
                  Create(middle_line),
                  run_time = 1.5)

        hexagon_size = 0.6
        horizontal_spacing = 1.5 * hexagon_size
        vertical_spacing = np.sqrt(3) * hexagon_size

        num_rows = 4
        num_cols = 4
        hexagons_left = VGroup()
        fine_mesh_lines = VGroup()
        degrees = VGroup()

        red_colors = [interpolate_color(RED_A, RED_E, alpha) for alpha in np.linspace(0, 1, 11)]

        for i in range(num_rows):
            for j in range(num_cols):
                center_x = j * horizontal_spacing
                center_y = i * vertical_spacing

                if j % 2 == 1:
                    center_y += vertical_spacing / 2

                center = [center_x, center_y, 0]
                hexagon = RegularPolygon(n=6, radius=hexagon_size, color = WHITE)
                hexagon.move_to(center)
                hexagons_left.add(hexagon)
                
                degree = random.randint(0, 10)
                degree_str = str(degree)
                degree_tex = Tex(degree_str, color = red_colors[degree])
                degree_tex.move_to(center)
                degrees.add(degree_tex)
                
                if random.random() < 0.5:
                    vertices = hexagon.get_vertices()
                    for k in range(len(vertices)):
                        fine_mesh_lines.add(Line(center, vertices[k], color=BLUE_C, stroke_width=2))
                        
                        if k > 0 and random.random() < 0.25:
                            fine_mesh_lines.add(Line(center, (vertices[k] + vertices[k-1])/2, color=BLUE_B, stroke_width=2))
                
        # We have to move hexagons_left and fine_mesh_lines
        # together otherwise they might get misplaced
        left_and_lines = VGroup()
        left_and_lines.add(hexagons_left)
        left_and_lines.add(fine_mesh_lines)
        left_and_lines.move_to([-x_position, 0, 0])

        hexagons_right = hexagons_left.copy()
        hexagons_right.move_to([x_position, 0, 0])
        degrees.move_to([x_position, 0, 0])

        self.play(LaggedStart(
            Create(hexagons_left, lag_ratio=0.15),
                  Create(hexagons_right, lag_ratio=0.15),
                  lag_ratio = 0.2,
                  run_time = 2
        ))

        self.play(LaggedStart(
            Create(fine_mesh_lines), lag_ratio = 0.3),
            Write(degrees),
            lag_ratio = 0.2,
            run_time = 7.5
        )
        
        all_objects = VGroup(*self.mobjects)
        self.play(all_objects.animate.shift(UP).fade(1))
        self.wait(2)

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

class scene4(ThreeDScene):
    def construct(self):
        # === Part 1: Equation Display (Initial 2D Camera View) ===
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES)
        self.camera.frame_height = 8.0 # Standard viewing height

        # 1. Title
        title = Tex("Reaction-Diffusion Systems", font_size=48)
        title.to_edge(UP, buff=0.7) # Position at the top
        self.play(Write(title), run_time=1.5)
        self.wait(0.3)

        # 2. Heat Equation with Non-linear Reaction Term
        # Using D for the diffusion coefficient
        equation = MathTex(
            r"\frac{\partial u}{\partial t} - D \Delta u + R(u) = f",
            font_size=42
        )
        
        # Initial and Boundary Conditions text
        ic_text = MathTex(r"u(\mathbf{x}, 0) = u_0(\mathbf{x})", r"\quad (\text{Initial Condition})", font_size=36)
        bc_text = MathTex(r"\text{+ Appropriate Boundary Conditions on } \partial\Omega", font_size=36)
        
        # Group equation and conditions, then position below title
        equation_group = VGroup(equation, ic_text, bc_text).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
        equation_group.next_to(title, DOWN, buff=0.6)
        
        self.play(Write(equation), run_time=2.0)
        self.wait(0.2)
        self.play(Write(ic_text), run_time=1.5)
        self.wait(0.2)
        self.play(Write(bc_text), run_time=1.5)
        self.wait(2.0) # Hold for viewing

        # 3. Disappearance of text elements
        text_elements_to_fade = VGroup(title, equation_group)
        self.play(FadeOut(text_elements_to_fade, shift=UP * 0.3), run_time=1.0)
        self.wait(0.5) # Pause before 3D part

        # === Part 2: 3D Time-Dependent Wave ===
        # Transition camera to a 3D perspective
        self.move_camera(
            phi=60 * DEGREES,       # Tilt
            theta=-75 * DEGREES,    # Rotate
            distance=16,            # Zoom out for a "big" grid view
            run_time=2.0
        )
        # self.camera.frame_height = 9 # Optionally adjust frame height for new view

        # Domain parameters for the square [-L, L] x [-L, L]
        L = 2.0 

        # 4. "Big" 3D Cartesian Grid
        axes = ThreeDAxes(
            x_range=[-L * 1.1, L * 1.1, L / 2],  # Extend ranges slightly beyond L
            y_range=[-L * 1.1, L * 1.1, L / 2],
            z_range=[-0.5, 1.5, 0.5],      # Adjust z_range to fit wave amplitude
            x_length=7, y_length=7, z_length=4, # Visual size of axes
            axis_config={"include_numbers": True, "font_size": 20, "include_tip": False},
        )
        axes_labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u(x,y,t)")
        
        # 5. Square domain on XY-plane (optional visual)
        domain_square = Polygon(
            axes.c2p(-L, -L, 0), axes.c2p(L, -L, 0),
            axes.c2p(L, L, 0), axes.c2p(-L, L, 0),
            color=BLUE_D, stroke_width=2.5, fill_opacity=0.1
        )
        
        self.play(Create(axes), Write(axes_labels), Create(domain_square), run_time=2.0)
        self.wait(0.5)

        # 6. Time-dependent wave solution
        time = ValueTracker(0) # This will track the current time 't' for the animation

        # Define the wave function u(x, y, t_normalized)
        # t_normalized will go from 0 to 1 during the animation.
        def moving_gaussian_wave(x, y, t_normalized):
            amplitude = 1.2
            width_param = 0.6  # Controls the "spread" of the Gaussian pulse
            
            # Wave travels from x = -L (left) to x = +L (right) as t_normalized goes 0 to 1
            # Start slightly off-screen and end slightly off-screen for smooth entry/exit
            start_x_center = -L - 2 * width_param 
            end_x_center = L + 2 * width_param
            current_x_center = interpolate(start_x_center, end_x_center, t_normalized)
            
            center_y = 0 # Wave centered along y=0 for this example
            
            # Gaussian pulse shape
            exponent = -(((x - current_x_center)**2 + (y - center_y)**2) / (2 * width_param**2))
            return amplitude * np.exp(exponent)

        # Create the Surface mobject. It will be updated based on `time.get_value()`.
        solution_surface = Surface(
            lambda u, v: axes.c2p(u, v, moving_gaussian_wave(u, v, time.get_value())),
            u_range=[-L, L],  # Corresponds to x-axis domain
            v_range=[-L, L],  # Corresponds to y-axis domain
            resolution=(56, 56), # (nu, nv) - samples for smoothness
            fill_opacity=0.75,
            checkerboard_colors=[TEAL_D, PURPLE_D], # Colors for the surface
            stroke_width=0.2,
            stroke_color=BLACK
        )

        # Add an updater to the surface so it redraws each frame based on the 'time' ValueTracker
        solution_surface.add_updater(
            lambda mob: mob.become( # .become() regenerates the mobject
                Surface(
                    lambda u, v: axes.c2p(u, v, moving_gaussian_wave(u, v, time.get_value())),
                    u_range=[-L, L], v_range=[-L, L],
                    resolution=(56, 56), # Keep consistent resolution
                    fill_opacity=0.75,
                    checkerboard_colors=[TEAL_D, PURPLE_D],
                    stroke_width=0.2, stroke_color=BLACK
                )
            )
        )
        
        self.play(Create(solution_surface), run_time=1.0) # Initial appearance of the wave at t=0
        
        # Animate the wave by changing the 'time' ValueTracker from 0 to 1
        animation_run_time = 7.0 # Duration of the wave movement animation
        self.play(
            time.animate.set_value(1), # Animate 'time' from its current value (0) to 1
            run_time=animation_run_time,
            rate_func=linear # Wave moves at a constant speed
        )
        self.wait(1.0) # Hold the final state of the wave

        # Cleanup
        all_3d_elements = VGroup(axes, axes_labels, domain_square, solution_surface)
        self.play(FadeOut(all_3d_elements), run_time=1.5)
        self.wait(0.5)

class scena5(Scene):
    def construct(self):
        invisible = Text("invisible").to_corner(UP)
        x_delta = config.frame_width / 6
        start_pos_l = [-x_delta, invisible.get_y(), 0]
        end_pos_l = [-x_delta, -invisible.get_y(), 0]
        start_pos_r = [x_delta, invisible.get_y(), 0]
        end_pos_r = [x_delta, -invisible.get_y(), 0]
        left_line = Line(start=start_pos_l, end=end_pos_l)
        right_line = Line(start=start_pos_r, end=end_pos_r)

        t1 = Text("Shape optimality", color=BLUE_C, font_size = 34).to_corner(UP)
        t2 = Text("Optimal degree", color=RED_C, font_size = 34)
        t3 = Text("Derefinement", color=GREEN_C, font_size = 34)

        t1.move_to([-2*x_delta, t1.get_y(), 0])
        t2.move_to([0, t1.get_y(), 0])
        t3.move_to([2*x_delta, t1.get_y(), 0])

        self.play(Create(left_line),
                  Create(right_line),
                  Write(t1),
                  Write(t2),
                  Write(t3),
                  run_time = 2.5)
        
        pentagon = RegularPolygon(n=5, radius=1.2, color = WHITE).move_to([-2*x_delta, 0, 0])
        baseline = Line([-1.5, pentagon.get_bottom()[1], 0], [1.5, pentagon.get_bottom()[1], 0], color=WHITE)

        big_square_size = 2.5
        bigsquare = Square(side_length=big_square_size, color=WHITE)
        bigsquare.move_to([2*x_delta, big_square_size/2 + pentagon.get_bottom()[1], 0])

        num_divisions = 4
        square_lines = VGroup()
        small_square_side = bigsquare.side_length / num_divisions

        for i in range(num_divisions-1):
            delta = (i + 1) * small_square_side
            square_line = Line([bigsquare.get_left()[0] + delta, bigsquare.get_top()[1], 0],
                               [bigsquare.get_left()[0] + delta, bigsquare.get_bottom()[1], 0],
                               color = WHITE)
            square_lines.add(square_line)
            square_line = Line([bigsquare.get_left()[0], bigsquare.get_bottom()[1] + delta, 0],
                               [bigsquare.get_right()[0], bigsquare.get_bottom()[1] + delta, 0],
                               color = WHITE)
            square_lines.add(square_line)

        self.play(Create(pentagon),
                  Create(baseline),
                  Create(bigsquare),
                  Create(square_lines),
                  run_time = 3)
        
        square_size = 2.4*math.sin(math.pi/5)
        square = Square(side_length=square_size, color=WHITE)
        square.move_to([-2*x_delta, square_size/2 + pentagon.get_bottom()[1], 0])

        square_vert = square.get_vertices()
        penta_vert = pentagon.get_vertices()

        # Vertices are ordered counter-clockwise from the top
        line1_1 = Line(square_vert[1], penta_vert[1], color=BLUE_A)
        line1_2 = Line(square_vert[1], penta_vert[0], color=BLUE_A)
        line1_3 = Line(square_vert[0], penta_vert[0], color=BLUE_A)
        line1_4 = Line(square_vert[0], penta_vert[4], color=BLUE_A)

        mesh_refinement1 = VGroup()
        mesh_refinement1.add(square)
        mesh_refinement1.add(line1_1)
        mesh_refinement1.add(line1_2)
        mesh_refinement1.add(line1_3)
        mesh_refinement1.add(line1_4)

        def poly1(x):
            return x
        def poly2(x):
            return x**2
        def poly3(x):
            return x**3
        
        axes = Axes(
            x_range=[-1.8, 1.8, 0.1], 
            y_range=[0, 5, 0.1],
            tips=False
        ).scale(0.25).move_to([0, pentagon.get_y(), 0] + UP/2)

        indexes = random.sample(range(0, len(square_lines)), 3)
        dashed_1 = DashedLine(square_lines[indexes[0]].get_start(), square_lines[indexes[0]].get_end(), dashed_ratio=0.5, color = GREEN_C)

        x_range = [-1.8, 1.8]
        self.play(Create(mesh_refinement1),
                  Create(axes.plot(poly1, x_range=x_range, color=RED_A)),
                  FadeOut(square_lines[indexes[0]]),
                  Create(dashed_1),
            run_time = 3)
        
        center = pentagon.get_center()
        line2_1 = Line(center, penta_vert[0], color=BLUE_B)
        line2_2 = Line(center, penta_vert[1], color=BLUE_B)
        line2_3 = Line(center, penta_vert[2], color=BLUE_B)
        line2_4 = Line(center, penta_vert[3], color=BLUE_B)
        line2_5 = Line(center, penta_vert[4], color=BLUE_B)

        mesh_refinement2 = VGroup()
        mesh_refinement2.add(line2_1)
        mesh_refinement2.add(line2_2)
        mesh_refinement2.add(line2_3)
        mesh_refinement2.add(line2_4)
        mesh_refinement2.add(line2_5)

        dashed_2 = DashedLine(square_lines[indexes[1]].get_start(), square_lines[indexes[1]].get_end(), dashed_ratio=0.5, color = GREEN_C)

        self.play(Transform(mesh_refinement1, mesh_refinement2),
                  Create(axes.plot(poly2, x_range=x_range, color=RED_B)),
                  FadeOut(square_lines[indexes[1]]),
                  Create(dashed_2),
            run_time = 3)
        
        pentagon2 = RegularPolygon(n=5, radius=0.6, color = BLUE_C).move_to([-2*x_delta, 0, 0])
        penta_vert2 = pentagon2.get_vertices()

        line3_1 = Line(penta_vert[0], penta_vert2[0], color=BLUE_C)
        line3_2 = Line(penta_vert[1], penta_vert2[1], color=BLUE_C)
        line3_3 = Line(penta_vert[2], penta_vert2[2], color=BLUE_C)
        line3_4 = Line(penta_vert[3], penta_vert2[3], color=BLUE_C)
        line3_5 = Line(penta_vert[4], penta_vert2[4], color=BLUE_C)

        mesh_refinement3 = VGroup()
        mesh_refinement3.add(pentagon2)
        mesh_refinement3.add(line3_1)
        mesh_refinement3.add(line3_2)
        mesh_refinement3.add(line3_3)
        mesh_refinement3.add(line3_4)
        mesh_refinement3.add(line3_5)

        dashed_3 = DashedLine(square_lines[indexes[2]].get_start(), square_lines[indexes[2]].get_end(), dashed_ratio=0.5, color = GREEN_C)

        self.play(Transform(mesh_refinement1, mesh_refinement3),
                  Create(axes.plot(poly3, x_range=x_range, color=RED_C)),
                  FadeOut(square_lines[indexes[2]]),
                  Create(dashed_3),
            run_time = 3)
        
        all_objects = VGroup(*self.mobjects)
        self.play(all_objects.animate.shift(UP).fade(1))
        self.wait(2)
