from manim import *
import numpy as np
import random
import math

class scene1(Scene):
    def construct(self):
        title_text = Tex("hp-adaptive strategies for high-fidelity multi-physics simulations",
                         font_size=40, color=WHITE)
        self.play(Write(title_text), run_time=2)
        self.wait(0.5)

        pietro = Text("Pietro Fumagalli", font_size=30, color=BLUE_C).shift(1.8 * LEFT + 0.5 * DOWN)
        francesco = Text("Francesco Derme", font_size=30, color=BLUE_C).shift(1.8 * RIGHT + 0.5 * DOWN)

        self.play(Write(pietro), 
                  Write(francesco),
                  run_time=1.5)
        
        self.wait(2)

        main_titles_group = VGroup(title_text, pietro, francesco)
        self.play(FadeOut(main_titles_group, shift=UP*0.5), run_time=1.0)
        self.wait(0.2)

        hexagons = VGroup()
        nrows = 3
        ncols = 6
        radius = 2
        startx = -config.frame_width/2
        starty = config.frame_height/2

        for i in range(nrows):
            for j in range(ncols):
                deltax = j * radius * 1.5
                deltay = -i * math.sqrt(3) * 2

                if j%2==1:
                    deltay -= radius * 0.5 * math.sqrt(3)

                hexagon = RegularPolygon(n=6, radius=radius, color=WHITE,
                                         stroke_width=4,
                                         stroke_opacity=random.random())
                hexagon.move_to([startx + deltax, starty + deltay, 0])
                hexagons.add(hexagon)

        self.play(LaggedStart(
            *[Create(hexagon) for hexagon in hexagons],
            lag_ratio=0.2),
            runt_time = 2.5)
        self.wait(0.5)

        coarse_mesh_lines = VGroup()
        indexes = [0, 5, 7, 8, 16]
        for ind in indexes:
            center = hexagons[ind].get_center()
            vertices = hexagons[ind].get_vertices()
            for vertex in vertices:
                coarse_mesh_lines.add(Line(center, vertex, color=BLUE_C, stroke_width=2.5))
        
        self.play(LaggedStart(
            *[Create(line) for line in coarse_mesh_lines],
            lag_ratio=0.2),
            run_time=2.0)
        self.wait(0.5)

        fine_mesh_lines = VGroup()
        indexes = [8, 16]
        for ind in indexes:
            center = hexagons[ind].get_center()
            vertices = hexagons[ind].get_vertices()
            for i in range(6):
                v1 = vertices[i]
                v2 = vertices[(i + 1) % 6]
                m01 = (center + v1) / 2
                m12 = (v1 + v2) / 2 
                m02 = (center + v2) / 2
                fine_mesh_lines.add(DashedLine(m01, m12, dashed_ratio=0.5, color = GREEN_C, stroke_width=1.5))
                fine_mesh_lines.add(DashedLine(m12, m02, dashed_ratio=0.5, color = GREEN_C, stroke_width=1.5))
                fine_mesh_lines.add(DashedLine(m02, m01, dashed_ratio=0.5, color = GREEN_C, stroke_width=1.5))
        
        self.play(LaggedStart(
            *[Create(line) for line in fine_mesh_lines],
            lag_ratio=0.05),
            run_time=2.0)
        self.wait(2)
        
        all_objects = VGroup(hexagons, coarse_mesh_lines, fine_mesh_lines)
        self.play(FadeOut(all_objects, shift=UP*0.5), run_time=1.0)
        self.wait(0.5)

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
            run_time = 2,
            lag_ratio = 0.3,
        ))

        self.play(LaggedStart(
            *[Create(line) for line in fine_mesh_lines],
            lag_ratio=0.2),
            Write(degrees),
            lag_ratio = 0.2,
            run_time = 6.5
        )
        
        self.wait(2.5)
        all_objects = VGroup(hexagons_left, hexagons_right, h, p, middle_line,
                             fine_mesh_lines, degrees)
        self.play(FadeOut(all_objects, shift=UP*0.5), run_time=1.0)
        self.wait(0.5)

class scene3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES)
        self.camera.frame_height = 8.5

        title = Text("A first example: the Laplacian", font_size=40).to_edge(UP).shift(DOWN * 0.5)
        self.play(Write(title), run_time=1.5)
        self.wait(0.5)

        laplace_eq_system = MathTex(
            r"""\begin{cases}
                -\Delta u(\mathbf{x}) = f(\mathbf{x}) & \mathbf{x} \in \Omega \\
                u(\mathbf{x}) = g(\mathbf{x}) & \mathbf{x} \in \partial\Omega 
            \end{cases}""",
            font_size=40
        ).next_to(title, DOWN, buff=0.35)

        self.play(Write(laplace_eq_system), run_time=2.5)

        boundary_expr = MathTex(
            r"""\text{Let } \Omega = (-L, L)^2, \partial\Omega = \bigcup_{i=1}^4 \Gamma_i \text{ where:}""",
            font_size=40
        ).next_to(laplace_eq_system, 3 * DOWN)

        self.play(Write(boundary_expr), run_time=2)
        
        boundary_details_1 = MathTex(            
            r"""\begin{alignedat}{2}
                \Gamma_1 &= \{ (L, y)   &&\mid -L \le y \le L \}, \text{ }
                \Gamma_2 &= \{ (-L, y)  &&\mid -L \le y \le L \}
            \end{alignedat}""",
            font_size=40
        ).next_to(boundary_expr, DOWN, buff = 0.5)

        boundary_details_2 = MathTex(            
            r"""\begin{alignedat}{2}
                \Gamma_3 &= \{ (x, L)   &&\mid -L \le x \le L \}, \text{ }
                \Gamma_4 &= \{ (x, -L)  &&\mid -L \le x \le L \}
            \end{alignedat}""",
            font_size=40
        ).next_to(boundary_details_1, DOWN)
     
        self.play(
            Write(boundary_details_1),
            Write(boundary_details_2),
            run_time=2)
        self.wait(2.0)

        text_elements = VGroup(title, laplace_eq_system, boundary_expr,
                               boundary_details_1, boundary_details_2)
        self.play(FadeOut(text_elements, shift=UP * 0.5), run_time=1.0)

        L_val = 1.0
        axes = ThreeDAxes(
            x_range=[-L_val*1.2, L_val*1.2, 0.2], 
            y_range=[-L_val*1.2, L_val*1.2, 0.2], 
            z_range=[-0.5, 1.5, 0.2],
            x_length=5.5, y_length=5.5, z_length=3.8,
            axis_config={"include_numbers": False, "font_size": 18, "include_tip": True},
        ).shift(1.5 * LEFT + 2.5 * DOWN)

        square = Square(side_length=4).shift(3.5 * RIGHT)
        half_side = square.side_length / 2

        vert_middle = Line([square.get_left()[0] + half_side, square.get_top()[1], 0],
                           [square.get_left()[0] + half_side, square.get_bottom()[1], 0],
                           color = WHITE)
        horiz_middle = Line([square.get_left()[0], square.get_top()[1] - half_side, 0],
                            [square.get_right()[0], square.get_top()[1] - half_side, 0],
                            color = WHITE)
        
        square_and_mid = VGroup(square, vert_middle, horiz_middle)
        self.add_fixed_in_frame_mobjects(square_and_mid)

        self.begin_ambient_camera_rotation(30*DEGREES, about='phi')
        self.begin_ambient_camera_rotation(20*DEGREES, about='theta')
        axes_labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u")
        self.play(LaggedStart(
            Create(axes),
            Write(axes_labels),
            Create(square_and_mid),
            lag_ratio = 0.2,
            run_time=2))
        self.stop_ambient_camera_rotation(about='phi')
        self.stop_ambient_camera_rotation(about='theta')

        divisions_per_quadrant = 4
        square_lines = VGroup()
        coarse_square_side = half_side / divisions_per_quadrant
        fine_square_side = coarse_square_side / 2

        colored_squares = VGroup()
        colors = [interpolate_color(RED_C, BLUE_C, alpha) for alpha in np.linspace(0, 1, 10)]

        def add_lines_and_squares(num, side, up_to_down, stroke_width):
            for i in range(num):
                delta = (i + 1) * side
                left = square.get_left()[0]
                top = square.get_top()[1]

                if up_to_down == -1:
                    top -= half_side

                if i != num-1:
                    # Vertical
                    square_line = Line([left + delta, top, 0], [left + delta, top - half_side, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)
                    # Horizontal
                    square_line = Line([left, top - delta, 0], [left + half_side, top - delta, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)
                    # Vertical shifted
                    square_line = Line([left + half_side + delta, top - up_to_down * half_side, 0],
                                    [left + half_side + delta, top - half_side - up_to_down * half_side, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)
                    # Horizontal shifted
                    square_line = Line([left + half_side, top - delta - up_to_down * half_side, 0],
                                    [left + 2 * half_side, top - delta - up_to_down * half_side, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)

                for j in range(num):
                    colored_square = Square(side_length=side)
                    colored_square.move_to([left + j * side + side/2, top - i * side - side/2, 0])
                    colored_square.set_stroke(width=0)
                    colored_square.set_fill(color=colors[random.randint(0, 9)], opacity=0.85)
                    colored_squares.add(colored_square)

                    colored_square = Square(side_length=side)
                    colored_square.move_to([left + j * side + side/2 + half_side,
                                            top - i * side - side/2 - up_to_down * half_side, 0])
                    colored_square.set_stroke(width=0)
                    colored_square.set_fill(color=colors[random.randint(0, 9)], opacity=0.85)
                    colored_squares.add(colored_square)

        add_lines_and_squares(divisions_per_quadrant, coarse_square_side, 1, 3)
        add_lines_and_squares(2*divisions_per_quadrant, fine_square_side, -1, 1.5)

        self.add_fixed_in_frame_mobjects(square_lines)
        self.add_fixed_in_frame_mobjects(colored_squares)

        def piecewise_solution_func(x, y):
            return math.sin(math.pi*x) * math.sin(math.pi*y) * math.exp(-7.5 * (x - y) * (x-y))
        
        surface_options = {
            "fill_opacity": 0.85,
            "checkerboard_colors": [BLUE_E, BLUE_C],
            "stroke_color": BLACK,
            "stroke_width": 0.1
        }

        solution_BL = Surface(
            lambda u, v: axes.c2p(u, v, piecewise_solution_func(u, v)),
            u_range=[-L_val, 0], v_range=[-L_val, 0],
            resolution=(18, 18),
            **surface_options
        )
        solution_BR = Surface(
            lambda u, v: axes.c2p(u, v, piecewise_solution_func(u, v)),
            u_range=[0, L_val], v_range=[-L_val, 0],
            resolution=(9, 9),
            **surface_options 
        )
        solution_TL = Surface(
            lambda u, v: axes.c2p(u, v, piecewise_solution_func(u, v)),
            u_range=[-L_val, 0], v_range=[0, L_val],
            resolution=(9, 9),
            **surface_options
        )
        solution_TR = Surface(
            lambda u, v: axes.c2p(u, v, piecewise_solution_func(u, v)),
            u_range=[0, L_val], v_range=[0, L_val],
            resolution=(18, 18),
            **surface_options 
        )
        
        solution = VGroup(solution_BL, solution_BR, solution_TL, solution_TR)

        self.play(
            Create(solution),
            Create(colored_squares),
            Create(square_lines),
            run_time=5.0)
        self.wait(4.0)

        all_3d_objects = VGroup(axes, axes_labels, solution)
        all_2d_objects = VGroup(square_and_mid, square_lines, colored_squares)
        self.play(FadeOut(all_3d_objects, shift=OUT*0.5),
                  FadeOut(all_2d_objects, shift=UP*0.5),
                  run_time=1.0)
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


class scene41(Scene): # Or scene41 if you prefer to keep your naming
    def construct(self):
        # --- 1. Domain and Detailed Non-homogeneous Mesh (Initially Centered) ---
        domain_size = 6.0
        L = domain_size / 2.0
        domain_outline = Square(side_length=domain_size, stroke_color=BLUE_D, stroke_width=2.5)
        domain_outline.move_to(ORIGIN)

        mesh_elements = VGroup()
        quadrant_L = L / 2.0
        quadrant_centers = {
            "TL": np.array([-quadrant_L, quadrant_L, 0]), "TR": np.array([quadrant_L, quadrant_L, 0]),
            "BL": np.array([-quadrant_L, -quadrant_L, 0]), "BR": np.array([quadrant_L, -quadrant_L, 0])
        }
        
        n_fine_divs = 12 
        n_coarse_divs = 6

        def create_subgrid_elements(center_pos, half_L_quad, num_divs, base_color):
            sub_elements = VGroup()
            step = (2 * half_L_quad) / num_divs
            for i in range(num_divs):
                for j in range(num_divs):
                    x_bottom_left = center_pos[0] - half_L_quad + i * step
                    y_bottom_left = center_pos[1] - half_L_quad + j * step
                    rect = Rectangle(width=step, height=step, 
                                     stroke_color=base_color, stroke_width=0.5,
                                     fill_opacity=0.0) 
                    rect.move_to(np.array([x_bottom_left + step/2, y_bottom_left + step/2, 0]))
                    sub_elements.add(rect)
            return sub_elements

        mesh_elements.add(*create_subgrid_elements(quadrant_centers["BL"], quadrant_L, n_fine_divs, DARK_GREY))
        mesh_elements.add(*create_subgrid_elements(quadrant_centers["TR"], quadrant_L, n_fine_divs, DARK_GREY))
        mesh_elements.add(*create_subgrid_elements(quadrant_centers["TL"], quadrant_L, n_coarse_divs, GREY_BROWN))
        mesh_elements.add(*create_subgrid_elements(quadrant_centers["BR"], quadrant_L, n_coarse_divs, GREY_BROWN))

        self.play(Create(domain_outline), Create(mesh_elements, lag_ratio=0.002, run_time=4.0))
        self.wait(0.5)
        mesh_display_group = VGroup(domain_outline, mesh_elements)

        # --- 2. Indicator Function I(x) Visualization on Mesh ---


        peak_bl_center = quadrant_centers["BL"]
        peak_tr_center = quadrant_centers["TR"]
        indicator_sigma = quadrant_L * 0.85

        def get_indicator_value(point_xyz):
            p = point_xyz[:2] 
            val_bl = 0.9 * np.exp(-np.sum((p - peak_bl_center[:2])**2) / (2 * indicator_sigma**2))
            val_tr = 0.9 * np.exp(-np.sum((p - peak_tr_center[:2])**2) / (2 * indicator_sigma**2))
            return np.clip(val_bl + val_tr, 0, 1.0) 

        indicator_value_animations = [
            element.animate.set_fill(
                interpolate_color(BLUE_E, RED_D, get_indicator_value(element.get_center())), 
                opacity=0.2 + 0.6 * get_indicator_value(element.get_center()) 
            ) for element in mesh_elements
        ]
        self.play(LaggedStart(*indicator_value_animations, lag_ratio=0.002), run_time=3.0)
        self.wait(1.0)

        # --- 3. Layout Change & 1D Plot of I(x) with Threshold ---
        target_mesh_scale = 0.75 
        target_mesh_position = LEFT * (self.camera.frame_width / 2) * 0.50 
        
        plot_width = self.camera.frame_width * 0.38 
        plot_height = 3.0
        plot_group_target_x = self.camera.frame_width / 3.8 
        plot_group_target_y = 0.2 

        plot_axes = Axes(
            x_range=[-L - 0.1, L + 0.1, L/2], 
            y_range=[-0.05, 1.05, 0.2],   
            x_length=plot_width,
            y_length=plot_height,
            axis_config={"include_numbers": True, "font_size": 16, 
                         "decimal_number_config": {"num_decimal_places": 1}},
            y_axis_config={"include_numbers": True, "font_size": 16,
                           "decimal_number_config": {"num_decimal_places": 1}},
            tips=False
        )
        plot_xlabel = MathTex("x \\text{ (domain slice)}", font_size=20).next_to(plot_axes.x_axis.get_right(), DR, buff=0.05)
        plot_ylabel = MathTex("I(x)", font_size=20).next_to(plot_axes.y_axis.get_top(), UL, buff=0.05)

        y_slice_for_plot = 0 
        indicator_curve = plot_axes.plot(
            lambda x_val: get_indicator_value(np.array([x_val, y_slice_for_plot, 0])),
            x_range=[-L, L],
            color=YELLOW_D, use_smoothing=True 
        )
        # --- CORRECTED LINE FOR CURVE_LABEL ---
        curve_label = Tex("$I(x,y_0)$", font_size=18, color=YELLOW_D).next_to(
            indicator_curve, 
            UR, # Direction is Upper-Right
            buff=0.15 
        )
        # --- END CORRECTION ---

        threshold_value_num = 0.35 
        threshold_line_on_plot = plot_axes.plot(lambda x_val: threshold_value_num, x_range=[-L,L], color=RED_B, stroke_width=3.5)
        threshold_label_on_plot = MathTex(f"T = {threshold_value_num}", font_size=20, color=RED_B)
        threshold_label_on_plot.next_to(threshold_line_on_plot.get_right(), RIGHT, buff=0.1)

        plot_group = VGroup(plot_axes, plot_xlabel, plot_ylabel, indicator_curve, curve_label, threshold_line_on_plot, threshold_label_on_plot)
        plot_group.move_to(np.array([plot_group_target_x, plot_group_target_y, 0]))

        self.play(
            mesh_display_group.animate.scale(target_mesh_scale).move_to(target_mesh_position),
            FadeIn(plot_group, shift=RIGHT * 0.2),
            run_time=1.5
        )
        self.wait(0.5)


        num_total_elements = len(mesh_elements.submobjects)
        example_indices = []
        if num_total_elements > 0:
            num_fine_elements_per_quad = n_fine_divs**2
            if num_fine_elements_per_quad > 0 : example_indices.append(num_fine_elements_per_quad // 2) 
            if num_total_elements > 2 * num_fine_elements_per_quad:
                 example_indices.append(2*num_fine_elements_per_quad + (n_coarse_divs**2 // 2) )
            if not example_indices and num_total_elements > 0: 
                 example_indices = [0, min(1, num_total_elements-1)] if num_total_elements > 1 else [0]

        CLUSTER_COLOR, NON_CLUSTER_COLOR = ORANGE, BLUE_C
        CLUSTER_OPACITY, NON_CLUSTER_OPACITY = 0.9, 0.25

        for element_idx in example_indices:
            if element_idx >= num_total_elements: continue 
            element = mesh_elements.submobjects[element_idx]
            value = get_indicator_value(element.get_center())
            
            original_stroke_color = element.get_stroke_color() 
            original_stroke_width = element.get_stroke_width()

            dot_on_plot = Dot(plot_axes.c2p(element.get_center()[0], value), color=YELLOW_A, radius=0.07)
            
            self.play(element.animate.set_stroke(YELLOW_A, width=4), Create(dot_on_plot), run_time=0.5)
            
            status_text_anims = []
            element_final_anim = None
            if value > threshold_value_num:
                status_text_obj = MathTex("I > T", color=GREEN_B, font_size=20).next_to(dot_on_plot, UR, buff=0.05)
                element_final_anim = element.animate.set_fill(CLUSTER_COLOR, opacity=CLUSTER_OPACITY).set_stroke(RED_E, width=1.5)
            else:
                status_text_obj = MathTex("I \\le T", color=RED_C, font_size=20).next_to(dot_on_plot, DR, buff=0.05)
                element_final_anim = element.animate.set_fill(NON_CLUSTER_COLOR, opacity=NON_CLUSTER_OPACITY).set_stroke(GREY, width=0.5)
            
            self.play(Write(status_text_obj), Indicate(threshold_line_on_plot, color=WHITE, scale_factor=1.05), run_time=0.8)
            if element_final_anim:
                self.play(element_final_anim, run_time=0.7)
            self.play(FadeOut(dot_on_plot), FadeOut(status_text_obj), 
                      element.animate.set_stroke(original_stroke_color, width=original_stroke_width),
                      run_time=0.5)
        
        final_element_animations = []
        processed_example_elements = [mesh_elements.submobjects[i] for i in example_indices if i < num_total_elements]

        for element in mesh_elements.submobjects:
            if element in processed_example_elements: 
                # Ensure its final state matches the last animation (if stroke was reset above)
                value = get_indicator_value(element.get_center())
                if value > threshold_value_num:
                     final_element_animations.append(element.animate.set_fill(CLUSTER_COLOR, opacity=CLUSTER_OPACITY).set_stroke(RED_E, width=1.5))
                else:
                     final_element_animations.append(element.animate.set_fill(NON_CLUSTER_COLOR, opacity=NON_CLUSTER_OPACITY).set_stroke(GREY, width=0.5))
                continue 

            value = get_indicator_value(element.get_center())
            if value > threshold_value_num:
                final_element_animations.append(element.animate.set_fill(CLUSTER_COLOR, opacity=CLUSTER_OPACITY).set_stroke(RED_E, width=1.0))
            else:
                final_element_animations.append(element.animate.set_fill(NON_CLUSTER_COLOR, opacity=NON_CLUSTER_OPACITY).set_stroke(GREY, width=0.5))
    
        
        self.wait(1.0)

        # --- Fade Out All ---
        self.play(*[FadeOut(mob) for mob in self.mobjects if mob is not None])
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
