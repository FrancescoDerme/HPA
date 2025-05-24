from manim import *
import numpy as np
import random
import math

class scene1(Scene):
    def construct(self):
        title_text = Tex("hp-adaptive strategies for high-fidelity multi-physics simulations",
                         font_size=40, color=WHITE)
        
        self.wait(1)
        self.play(Write(title_text), run_time=3)
        self.wait(1.5)

        pietro = Text("Pietro Fumagalli", font_size=30, color=BLUE_C).shift(1.8 * LEFT + 0.5 * DOWN)
        francesco = Text("Francesco Derme", font_size=30, color=BLUE_C).shift(1.8 * RIGHT + 0.5 * DOWN)

        self.play(Write(pietro), 
                  Write(francesco),
                  run_time=1.5)
        
        self.wait(5)

        main_titles_group = VGroup(title_text, pietro, francesco)
        self.play(FadeOut(main_titles_group, shift=UP*0.5), run_time=1.0)
        self.wait(0.5)

        """
        # Removed complicated transition
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
        """
        
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

        self.wait(1)
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
        
        self.wait(4.5)
        all_objects = VGroup(hexagons_left, hexagons_right, h, p, middle_line,
                             fine_mesh_lines, degrees)
        self.play(FadeOut(all_objects, shift=UP*0.5), run_time=1.0)
        self.wait(0.5)

class scene3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES)
        self.camera.frame_height = 8.5

        self.wait(1)
        title = Text("A first example: the laplacian", font_size=40).to_edge(UP).shift(DOWN * 0.5)
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
        axes_labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u")

        square = Square(side_length=4, z_index = 1).shift(3.5 * RIGHT)
        half_side = square.side_length / 2
        left = square.get_left()[0]
        right = square.get_right()[0]
        top = square.get_top()[1]
        bot = square.get_bottom()[1]

        # z_index is needed to consistenly render
        # the lines on top of the colored squares.
        # This DOES NOT change the z coordinate.
        vert_middle = Line([left + half_side, top, 0], [left + half_side, bot, 0],
                           color = WHITE, z_index = 1)
        horiz_middle = Line([left, top - half_side, 0], [right, top - half_side, 0],
                            color = WHITE, z_index = 1)
        
        square_and_mid = VGroup(square, vert_middle, horiz_middle)
        self.add_fixed_in_frame_mobjects(square_and_mid)

        phi_tracker, theta_tracker, *_, = self.camera.get_value_trackers()

        self.play(
            phi_tracker.animate.set_value(60*DEGREES),
            theta_tracker.animate.set_value(-50*DEGREES),
            LaggedStart(
                Create(axes),
                Write(axes_labels),
                Create(square_and_mid),
                lag_ratio = 0.2,
            ),
            run_time=2
        )

        square_lines = VGroup()
        colored_squares = VGroup()

        divisions_per_quadrant = 4
        coarse_square_side = half_side / divisions_per_quadrant
        fine_square_side = coarse_square_side / 2
        colors = [interpolate_color(RED_C, BLUE_C, alpha) for alpha in np.linspace(0, 1, 10)]

        def add_lines_and_squares(num, side, up_to_down, stroke_width):
            for i in range(num):
                delta = (i + 1) * side
                quadrant_top = square.get_top()[1]

                if up_to_down == -1:
                    quadrant_top -= half_side

                if i != num - 1:
                    # Vertical
                    square_line = Line([left + delta, quadrant_top, 0], [left + delta, quadrant_top - half_side, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)
                    # Horizontal
                    square_line = Line([left, quadrant_top - delta, 0], [left + half_side, quadrant_top - delta, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)
                    # Vertical shifted
                    square_line = Line([left + half_side + delta, quadrant_top - up_to_down * half_side, 0],
                                    [left + half_side + delta, quadrant_top - half_side - up_to_down * half_side, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)
                    # Horizontal shifted
                    square_line = Line([left + half_side, quadrant_top - delta - up_to_down * half_side, 0],
                                    [left + 2 * half_side, quadrant_top - delta - up_to_down * half_side, 0],
                                    color = WHITE, stroke_width = stroke_width)
                    square_lines.add(square_line)

                for j in range(num):
                    colored_square = Square(side_length=side)
                    colored_square.move_to([left + j * side + side/2, quadrant_top - i * side - side/2, 0])
                    colored_square.set_stroke(width=0)
                    colored_square.set_fill(color=colors[random.randint(0, 9)], opacity=0.85)
                    colored_squares.add(colored_square)

                    colored_square = Square(side_length=side)
                    colored_square.move_to([left + j * side + side/2 + half_side,
                                            quadrant_top - i * side - side/2 - up_to_down * half_side, 0])
                    colored_square.set_stroke(width=0)
                    colored_square.set_fill(color=colors[random.randint(0, 9)], opacity=0.85)
                    colored_squares.add(colored_square)

        add_lines_and_squares(divisions_per_quadrant, coarse_square_side, 1, 3)
        add_lines_and_squares(2*divisions_per_quadrant, fine_square_side, -1, 1.5)

        # Order matters: we want to show lines on top of squares
        self.add_fixed_in_frame_mobjects(colored_squares)
        self.add_fixed_in_frame_mobjects(square_lines)

        def solution_func(x, y):
            return math.sin(math.pi*x) * math.sin(math.pi*y) * math.exp(-7.5 * (x - y) * (x-y))
        
        surface_options = {
            "fill_opacity": 0.85,
            "checkerboard_colors": [BLUE_E, BLUE_C],
            "stroke_color": BLACK,
            "stroke_width": 0.1
        }

        solution_BL = Surface(
            lambda u, v: axes.c2p(u, v, solution_func(u, v)),
            u_range=[-L_val, 0],
            v_range=[-L_val, 0],
            resolution=(24, 24),
            **surface_options
        )
        solution_BR = Surface(
            lambda u, v: axes.c2p(u, v, solution_func(u, v)),
            u_range=[0, L_val],
            v_range=[-L_val, 0],
            resolution=(9, 9),
            **surface_options 
        )
        solution_TL = Surface(
            lambda u, v: axes.c2p(u, v, solution_func(u, v)),
            u_range=[-L_val, 0],
            v_range=[0, L_val],
            resolution=(9, 9),
            **surface_options
        )
        solution_TR = Surface(
            lambda u, v: axes.c2p(u, v, solution_func(u, v)),
            u_range=[0, L_val],
            v_range=[0, L_val],
            resolution=(24, 24),
            **surface_options 
        )
        
        solution = VGroup(solution_BL, solution_BR, solution_TL, solution_TR)
        solution.set_z_index(1)

        self.play(
            Create(solution),
            Create(colored_squares),
            Create(square_lines),
            run_time=6.0)
        self.wait(4.0)

        all_3d_objects = VGroup(axes, axes_labels, solution)

        # Order matters: we want to fade square_and_mid before colored_squares
        all_2d_objects = VGroup(square_and_mid, colored_squares, square_lines)
        self.play(FadeOut(all_3d_objects, shift=OUT*0.5),
                  FadeOut(all_2d_objects, shift=UP*0.5),
                  run_time=1.0)
        self.wait(0.5)

class scene4(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES)
        self.camera.frame_height = 8.5

        title = Text("Time-dependent reaction-diffusion systems", font_size=40).to_edge(UP).shift(DOWN * 0.5)

        self.wait(1)
        self.play(Write(title), run_time=1.5)
        self.wait(0.5)

        rd_eq_system = MathTex(
            r"""\begin{cases}
                \frac{\partial u}{\partial t} - \mu \Delta u + \sigma u = f & (\mathbf{x}, \mathbf{t}) \in \Omega \times (0, T) \\
                u(\mathbf{x}, 0) = u_0(\mathbf{x}) & \mathbf{x} \in \Omega \\
                u(\mathbf{x}, \mathbf{t}) = g(\mathbf{x}, \mathbf{t}) & (\mathbf{x}, \mathbf{t}) \in \partial\Omega \times (0, T)
            \end{cases}""",
            font_size=40
        ).next_to(title, DOWN, buff=0.35)
                        
        self.play(Write(rd_eq_system), run_time=3.5)
        self.wait(3.0)

        text_elements = VGroup(title, rd_eq_system)
        self.play(FadeOut(text_elements, shift=UP * 0.5), run_time=1.0)
        self.wait(0.5)

        L_val = 1.0 
        axes = ThreeDAxes(
            x_range=[-L_val*1.2, L_val*1.2, 0.2],
            y_range=[-L_val*1.2, L_val*1.2, 0.2],
            z_range=[-0.5, 1.5, 0.2],
            x_length=5.5, y_length=5.5, z_length=3.8,
            axis_config={"include_numbers": False, "font_size": 18, "include_tip": True},
        ).shift(1.5 * LEFT + 2.5 * DOWN)
        axes_labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u")

        square = Square(side_length=4, z_index = 1).shift(3.5 * RIGHT)
        half_side = square.side_length / 2
        square_center_x = square.get_center()[0]
        square_center_y = square.get_center()[1]
        left = square.get_left()[0]
        right = square.get_right()[0]
        top = square.get_top()[1]
        bot = square.get_bottom()[1]

        vert_middle = Line([left + half_side, top, 0], [left + half_side, bot, 0],
                           color = WHITE, z_index = 1)
        horiz_middle = Line([left, top - half_side, 0], [right, top - half_side, 0],
                            color = WHITE, z_index = 1)
        
        square_and_mid = VGroup(square, vert_middle, horiz_middle)
        self.add_fixed_in_frame_mobjects(square_and_mid)

        phi_tracker, theta_tracker, *_, = self.camera.get_value_trackers()
        
        self.play(
            phi_tracker.animate.set_value(60*DEGREES),
            theta_tracker.animate.set_value(-50*DEGREES),
            LaggedStart(
                Create(axes),
                Write(axes_labels),
                Create(square_and_mid),
                lag_ratio = 0.2,
            ),
            run_time=2
        )
        
        time = ValueTracker(-2.5)

        def solution_func(x, y, t):            
            return math.sin(math.pi*(x-t)) * math.sin(math.pi*y) * math.exp(-7.5 * (x - t - y) * (x - t -y))

        square_lines = VGroup()
        colored_squares = VGroup()

        divisions = 8
        coarse_square_side = square.side_length / divisions
        colors = [interpolate_color(RED_C, BLUE_C, alpha) for alpha in np.linspace(0, 1, 10)]
        def add_lines_and_squares(num, side, stroke_width):
            for i in range(num):
                delta = (i + 1) * side

                if i != num - 1:
                    # Vertical
                    square_line = Line([left + delta, top, 0], [left + delta, bot, 0],
                                    color = WHITE, stroke_width = stroke_width, z_index = 1)
                    square_lines.add(square_line)
                    # Horizontal
                    square_line = Line([left, top - delta, 0], [right, top - delta, 0],
                                    color = WHITE, stroke_width = stroke_width, z_index = 1)
                    square_lines.add(square_line)

                for j in range(num):
                    colored_square = Square(side_length=side)
                    colored_square.move_to([left + j * side + side/2, top - i * side - side/2, 0])
                    colored_square.set_stroke(width=0)
                    colored_square.set_fill(color=colors[random.randint(0, 9)], opacity=0.85)
                    colored_squares.add(colored_square)

        add_lines_and_squares(divisions, coarse_square_side, 3)

        adaptive_grid = VGroup()
        highly_refined = {}
        refined = {}

        def update_adaptive_grid(mob):
            t = time.get_value()
            new_side = coarse_square_side / 2
            new_grid_mobjects = VGroup()

            for i in range(divisions):
                for j in range(divisions):
                    current_center_x = square_center_x - half_side + (j + 0.5) * coarse_square_side
                    current_center_y = square_center_y - half_side + (i + 0.5) * coarse_square_side
                    
                    deltax = current_center_x - square_center_x
                    deltay = current_center_y - square_center_y
                    mapped_x = deltax * L_val / half_side
                    mapped_y = deltay * L_val / half_side

                    solution_height = solution_func(mapped_x, mapped_y, t)
                    used_colors = []  

                    if solution_height > 0.6:
                        if (i*divisions+j) in refined:
                            refined.pop(i*divisions + j)

                        check_for_colors = 0
                        if (i*divisions+j) in highly_refined:
                            check_for_colors = 1

                        colored_squares[i*divisions + j].set_fill(opacity=0)

                        vertices = [
                            [current_center_x - new_side, current_center_y, 0],
                            [current_center_x - new_side, current_center_y + new_side, 0],
                            [current_center_x, current_center_y + new_side, 0],
                            [current_center_x + new_side, current_center_y + new_side, 0],
                            [current_center_x + new_side, current_center_y, 0],
                            [current_center_x + new_side, current_center_y - new_side, 0],
                            [current_center_x, current_center_y - new_side, 0],
                            [current_center_x - new_side, current_center_y - new_side, 0]
                        ]
   
                        for k in range(8):
                            p = Polygon(vertices[k], [current_center_x,  current_center_y, 0], vertices[(k+1)%8])
                            p.set_stroke(width=0)

                            if check_for_colors:
                                p.set_fill(color=colors[highly_refined[i*divisions+j][k]], opacity=0.9)
                            else:
                                color_index = random.randint(0, 9)
                                p.set_fill(color=colors[color_index], opacity=0.9)
                                used_colors.append(color_index)

                            new_grid_mobjects.add(p)

                            line = Line([current_center_x,  current_center_y, 0],
                                        vertices[k], color = WHITE, stroke_width = 3)
                            new_grid_mobjects.add(line)

                        if check_for_colors == 0:
                            highly_refined[i*divisions+j] = used_colors
                    elif solution_height > 0.2:
                        if (i*divisions+j) in highly_refined:
                            highly_refined.pop(i*divisions + j)

                        check_for_colors = 0
                        if (i*divisions+j) in refined:
                            check_for_colors = 1

                        colored_squares[i*divisions + j].set_fill(opacity=0)

                        centers = [
                            [current_center_x - new_side/2, current_center_y - new_side/2, 0],
                            [current_center_x - new_side/2, current_center_y + new_side/2, 0],
                            [current_center_x + new_side/2, current_center_y - new_side/2, 0],
                            [current_center_x + new_side/2, current_center_y + new_side/2, 0]
                        ]

                        for k in range(4):
                            s = Square(side_length=new_side)
                            s.move_to(centers[k])
                            s.set_stroke(width=0)

                            if check_for_colors:
                                s.set_fill(color=colors[refined[i*divisions+j][k]], opacity=0.9)
                            else:
                                color_index = random.randint(0, 9)
                                s.set_fill(color=colors[color_index], opacity=0.9)
                                used_colors.append(color_index)

                            new_grid_mobjects.add(s)
                        
                        line = Line([current_center_x, current_center_y + new_side, 0],
                                    [current_center_x, current_center_y - new_side, 0],
                                    color = WHITE, stroke_width = 3)
                        new_grid_mobjects.add(line)
                        line = Line([current_center_x + new_side, current_center_y, 0],
                                    [current_center_x - new_side, current_center_y, 0],
                                    color = WHITE, stroke_width = 3)
                        new_grid_mobjects.add(line)

                        if check_for_colors == 0:
                            refined[i*divisions+j] = used_colors
                    else:
                        if (i*divisions+j) in highly_refined:
                            highly_refined.pop(i*divisions + j)
                        if (i*divisions+j) in refined:
                            refined.pop(i*divisions + j)
                        colored_squares[i*divisions+j].set_fill(opacity=0.85)
                        
            mob.become(new_grid_mobjects)
            self.add_fixed_in_frame_mobjects(mob)

        self.add_fixed_in_frame_mobjects(adaptive_grid)
        self.add_fixed_in_frame_mobjects(colored_squares)
        self.add_fixed_in_frame_mobjects(square_lines)

        adaptive_grid.add_updater(update_adaptive_grid)
        update_adaptive_grid(adaptive_grid)

        surface_options = {
            "u_range": [-L_val, L_val],
            "v_range": [-L_val, L_val],
            "resolution": (24, 24),
            "fill_opacity": 0.85,
            "checkerboard_colors": [BLUE_E, BLUE_C],
            "stroke_color": BLACK,
            "stroke_width": 0.1,
        }

        solution_surface = Surface(
            lambda u, v: axes.c2p(u, v, solution_func(u, v, time.get_value())),
            **surface_options
        )

        solution_surface.add_updater(
            lambda mob: mob.become(
                Surface(
                    lambda u, v: axes.c2p(u, v, solution_func(u, v, time.get_value())),
                    **surface_options
                )
            )
        )
        
        self.play(Create(solution_surface),
                  Create(colored_squares),
                  Create(square_lines),
                  run_time=4.0)
        
        self.play(
            time.animate.set_value(2.5),
            run_time=11.0,
            rate_func=linear
        )

        self.wait(1.0)

        # If we don't clear the updater it will try to bring
        # the surface's opacity back to 0.85
        solution_surface.clear_updaters()
        adaptive_grid.clear_updaters()

        all_3d_objects = VGroup(solution_surface, axes, axes_labels)
        all_2d_objects = VGroup(square_and_mid, square_lines, colored_squares, adaptive_grid)
        self.play(FadeOut(all_3d_objects, shift=OUT*0.5),
                  FadeOut(all_2d_objects, shift=UP*0.5),
                  run_time=1.0)
        self.wait(0.5)
        
class scene5(Scene):
    def construct(self):
        domain_size = 6.0
        L_val = domain_size / 2.0
        domain = Square(side_length=domain_size, stroke_color=WHITE, stroke_width=3)
        domain.move_to(ORIGIN)

        quadrant_centers = {
            "TL": np.array([-L_val/2, L_val/2, 0]),
            "TR": np.array([L_val/2, L_val/2, 0]),
            "BL": np.array([-L_val/2, -L_val/2, 0]),
            "BR": np.array([L_val/2, -L_val/2, 0])
        }
        
        n_fine_divs = 12 
        n_coarse_divs = 6

        mesh_elements = VGroup()

        def add_rectangles(center_pos, num_divs, color):
            step = L_val / num_divs
            for i in range(num_divs):
                for j in range(num_divs):
                    x_bottom_left = center_pos[0] - L_val/2 + i * step
                    y_bottom_left = center_pos[1] - L_val/2 + j * step
                    r = Rectangle(width=step, height=step, 
                                  stroke_color=color, stroke_width=0.5,
                                  fill_opacity=0.0) 
                    r.move_to([x_bottom_left + step/2, y_bottom_left + step/2, 0])
                    mesh_elements.add(r)

        add_rectangles(quadrant_centers["BL"], n_fine_divs, WHITE)
        add_rectangles(quadrant_centers["TR"], n_fine_divs, WHITE)
        add_rectangles(quadrant_centers["TL"], n_coarse_divs, WHITE)
        add_rectangles(quadrant_centers["BR"], n_coarse_divs, WHITE)

        start_time = -4
        time = ValueTracker(start_time)

        def solution_func(x, y, t):            
            return math.sin(math.pi*(x-t)) * math.sin(math.pi*y) * math.exp(-7.5 * (x - t - y) * (x - t -y))
        
        def update_mesh_element(mob):
            t = time.get_value()
            value = solution_func(mob.get_center()[0] / 6, mob.get_center()[1] / 6, time.get_value())
            value = (value + 0.11) / 1.11
            r = Rectangle(width=mob.width, height=mob.height, 
                        stroke_color=mob.get_stroke_color(),
                        stroke_width=mob.get_stroke_width(),
                        fill_opacity = min(2*(t-start_time), 1),
                        fill_color = interpolate_color(BLUE_C, RED_C, value))
            r.move_to(mob.get_center())
            mob.become(r)
                
        for element in mesh_elements:
            element.add_updater(update_mesh_element)
            update_mesh_element(element)

        self.wait(1)
        self.play(Create(domain),
                  Create(mesh_elements, lag_ratio=0.002, run_time=4.0))
        self.wait(0.5)
        
        mesh = VGroup(domain, mesh_elements)
        target_mesh_scale = 0.75 
        target_mesh_position = LEFT * self.camera.frame_width / 4 
        target_plot_position = RIGHT * self.camera.frame_width / 5

        plot_axes = Axes(
            x_range=[-1.1, 1.1, 0.2], 
            y_range=[-0.075, 1.1, 0.2],   
            x_length=2.5,
            y_length=1.5,
            axis_config={
                "include_numbers": False,
                "include_tip": True,
                "tick_size": 0.075,
                "tip_width": 0.15,
                "tip_height": 0.2
            },
        ).scale(2.5)

        plot_axes.move_to(target_plot_position)

        x_param = 0.0
        indicator_params = {
            "x_range": [-1, 1],
            "color": BLUE_C,
            "use_smoothing": False 
        }

        indicator_curve = plot_axes.plot(
            lambda y_val: solution_func(x_param, y_val, start_time),
            **indicator_params
        )

        plot_label = MathTex("I(0, y)", font_size=42, color = BLUE_C)
        plot_label.move_to([plot_axes.x_axis.get_left()[0] - plot_label.get_left()[0], plot_axes.y_axis.get_bottom()[1], 0])
        plot_label.shift(UP*0.45)

        threshold = 0.5
        threshold_line = plot_axes.plot(lambda x_val: threshold, x_range=[-1.1,1.1], color=RED_C, stroke_width=3.5)
        threshold_label = MathTex("threshold", font_size=42, color=RED_C)
        threshold_label.move_to(threshold_line.get_left())
        threshold_label.shift(UP*0.2).shift(RIGHT) 

        lines_group = VGroup(plot_label,
                             threshold_line, threshold_label)
        
        wait_to_plot = 2.55
        
        def update_indicator(mob):
            t = time.get_value() + 0.5
            p = plot_axes.plot(lambda y_val: solution_func(x_param, y_val, t),
                               **indicator_params)
            
            if t < start_time + wait_to_plot:
                p.set_opacity(0)
            else:
                p.set_opacity(min(0.35, 1.5*(t - (wait_to_plot + start_time))))
            mob.become(p)
        
        self.add(indicator_curve)
        indicator_curve.add_updater(update_indicator)

        time_value_animation = time.animate(
            run_time=18, 
            rate_func=linear
        ).set_value(2)

        mesh_sequence = Succession(
            Wait(1.5),
            mesh.animate(run_time = 2,).scale(target_mesh_scale).move_to(target_mesh_position)
        )
        
        plot_sequence = Succession(
            Wait(wait_to_plot),
            FadeIn(plot_axes, shift=UP*0.5, run_time = 1),
            Write(lines_group, run_time = 1.5),
        )

        self.play(
            time_value_animation,
            mesh_sequence,
            plot_sequence,
            rate_func=linear
        )

        self.wait(1.0)

        indicator_curve.clear_updaters()
        for element in mesh_elements:
            element.clear_updaters()

        all_2d_objects = VGroup(plot_axes, lines_group, mesh, indicator_curve)
        self.play(FadeOut(all_2d_objects, shift=UP*0.5),
                  run_time=1.0)
        self.wait(0.5)

class scene6(Scene):
    def construct(self):
        invisible = Text("invisible").to_corner(UP)
        x_delta = config.frame_width / 6
        start_pos_l = [-x_delta, invisible.get_y(), 0]
        end_pos_l = [-x_delta, -invisible.get_y(), 0]
        start_pos_r = [x_delta, invisible.get_y(), 0]
        end_pos_r = [x_delta, -invisible.get_y(), 0]
        left_line = Line(start=start_pos_l, end=end_pos_l)
        right_line = Line(start=start_pos_r, end=end_pos_r)

        t1 = Text("Shape optimality", color=BLUE_C, font_size = 34).to_corner(UP).shift(0.25*DOWN)
        t1.move_to([-2*x_delta, t1.get_y(), 0])

        t2 = Text("Optimal degree", color=RED_C, font_size = 34)
        t2.move_to([0, t1.get_y(), 0])

        t3 = Text("Derefinement", color=GREEN_C, font_size = 34)
        t3.move_to([2*x_delta, t1.get_y(), 0])
        
        # Group 1
        pentagon = RegularPolygon(n=5, radius=1.2, color = WHITE).move_to([-2*x_delta, 0, 0])

        self.wait(1)
        self.play(LaggedStart(
                Write(t1),
                Create(pentagon),
                lag_ratio=0.2),
            run_time = 2
        )
        
        square_size = 2.4*math.sin(math.pi/5)
        square = Square(side_length=square_size, color=WHITE)
        square.move_to([-2*x_delta, square_size/2 + pentagon.get_bottom()[1], 0])
        square_vert = square.get_vertices()
        penta_vert = pentagon.get_vertices()

        # Vertices are ordered counter-clockwise from the top
        mesh_refinement1 = VGroup()
        mesh_refinement1.add(square)
        mesh_refinement1.add(Line(square_vert[1], penta_vert[1], color=BLUE_A))
        mesh_refinement1.add(Line(square_vert[1], penta_vert[0], color=BLUE_A))
        mesh_refinement1.add(Line(square_vert[0], penta_vert[0], color=BLUE_A))
        mesh_refinement1.add(Line(square_vert[0], penta_vert[4], color=BLUE_A))

        self.play(Create(mesh_refinement1),
                  run_time = 2)
        
        mesh_refinement2 = VGroup()
        center = pentagon.get_center()
        mesh_refinement2.add(Line(center, penta_vert[0], color=BLUE_B))
        mesh_refinement2.add(Line(center, penta_vert[1], color=BLUE_B))
        mesh_refinement2.add(Line(center, penta_vert[2], color=BLUE_B))
        mesh_refinement2.add(Line(center, penta_vert[3], color=BLUE_B))
        mesh_refinement2.add(Line(center, penta_vert[4], color=BLUE_B))

        self.play(Transform(mesh_refinement1, mesh_refinement2),
                  run_time = 2)
        
        pentagon2 = RegularPolygon(n=5, radius=0.6, color = BLUE_C).move_to([-2*x_delta, 0, 0])
        penta_vert2 = pentagon2.get_vertices()

        mesh_refinement3 = VGroup()
        mesh_refinement3.add(Line(penta_vert[0], penta_vert2[0], color=BLUE_C))
        mesh_refinement3.add(Line(penta_vert[1], penta_vert2[1], color=BLUE_C))
        mesh_refinement3.add(Line(penta_vert[2], penta_vert2[2], color=BLUE_C))
        mesh_refinement3.add(Line(penta_vert[3], penta_vert2[3], color=BLUE_C))
        mesh_refinement3.add(Line(penta_vert[4], penta_vert2[4], color=BLUE_C))
        mesh_refinement3.add(pentagon2)

        self.play(Transform(mesh_refinement1, mesh_refinement3),
                  run_time = 2)

        # Group 2
        def poly1(x):
            return x
        def poly2(x):
            return x**2
        def poly3(x):
            return x**3
        
        axes = Axes(
            x_range=[-1, 1, 0.2], 
            y_range=[-0.2, 1, 0.2],
            axis_config={
                "include_numbers": False,
                "include_tip": True,
                "tick_size": 0.4,
                "tip_width": 0.6,
                "tip_height": 0.8
            },
        ).scale(0.275)
        
        axes.move_to([0, pentagon.get_bottom()[1] + axes.get_top()[1] - 0.2, 0])
        
        self.play(LaggedStart(
                Create(left_line),
                Write(t2),
                Create(axes),
                lag_ratio = 0.2),
            run_time = 2
        )
        
        x_range=[-1, 1]
        plot1 = axes.plot(poly1, x_range=x_range, color=RED_A)
        
        self.play(Create(plot1),
                  run_time = 2)
        
        plot2 = axes.plot(poly2, x_range=x_range, color=RED_B)

        self.play(Create(plot2),
                  run_time = 2)
        
        plot3 = axes.plot(poly3, x_range=x_range, color=RED_C)

        self.play(Create(plot3),
                  run_time = 2)
        
        # Group 3        
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

        self.play(LaggedStart(
                Create(right_line),
                Write(t3),
                Create(bigsquare),
                Create(square_lines),
                lag_ratio=0.2),
            run_time = 2
            )

        indexes = [2, 4, 5]
        dashed_1 = DashedLine(square_lines[indexes[0]].get_start(), square_lines[indexes[0]].get_end(), dashed_ratio=0.5, color = GREEN_C)

        self.play(FadeOut(square_lines[indexes[0]]),
                  Create(dashed_1),
                  run_time = 2)

        dashed_2 = DashedLine(square_lines[indexes[1]].get_start(), square_lines[indexes[1]].get_end(), dashed_ratio=0.5, color = GREEN_C)

        self.play(FadeOut(square_lines[indexes[1]]),
                  Create(dashed_2),
                  run_time = 2)

        dashed_3 = DashedLine(square_lines[indexes[2]].get_start(), square_lines[indexes[2]].get_end(), dashed_ratio=0.5, color = GREEN_C)

        self.play(FadeOut(square_lines[indexes[2]]),
                  Create(dashed_3),
                  run_time = 2)
        
        self.wait(2.5)
        
        square_lines.remove(square_lines[indexes[2]])
        square_lines.remove(square_lines[indexes[1]])
        square_lines.remove(square_lines[indexes[0]])

        all_2d_objects = VGroup(left_line, right_line, t1, t2, t3,
                                pentagon, bigsquare, square_lines,
                                dashed_1, dashed_2,
                                mesh_refinement1, dashed_3,
                                axes, plot1, plot2, plot3)
        self.play(FadeOut(all_2d_objects, shift = UP*0.5),
                  run_time=1.0)
        self.wait(0.5)

def _create_trophy_geometry() -> VGroup:
    cup_width_top = 0.8
    cup_width_bottom = 0.4
    cup_height = 0.6
    stem_width = 0.15
    stem_height = 0.5
    base_width = 0.7
    base_height = 0.12

    cup_points = [
        (-cup_width_top / 2, cup_height / 2, 0),
        (cup_width_top / 2, cup_height / 2, 0),
        (cup_width_bottom / 2, -cup_height / 2, 0),
        (-cup_width_bottom / 2, -cup_height / 2, 0)
    ]

    cup = Polygon(*cup_points)
    stem = Rectangle(width=stem_width, height=stem_height)
    base = Rectangle(width=base_width, height=base_height)

    stem.next_to(cup, DOWN, buff=0)
    base.next_to(stem, DOWN, buff=0)

    trophy_geometry = VGroup(cup, stem, base)
    trophy_geometry.move_to(ORIGIN, aligned_edge=DOWN)
    return trophy_geometry

def create_trophy_outline(outline_color=WHITE, stroke_width_val=4) -> VGroup:
    trophy_outline = _create_trophy_geometry()
    for part in trophy_outline:
        part.set_fill(opacity=0.0)
        part.set_stroke(color=outline_color, width=stroke_width_val)
    return trophy_outline

def create_filled_trophy_shape(fill_color=GOLD, stroke_color=None, stroke_width_val=0) -> VGroup:
    trophy_filled = _create_trophy_geometry()
    actual_stroke_color = stroke_color if stroke_color is not None else fill_color
    for part in trophy_filled:
        part.set_fill(color=fill_color, opacity=1.0)
        part.set_stroke(color=actual_stroke_color, width=stroke_width_val)
    return trophy_filled

class scene7(Scene):
    def construct(self):
        num_trophies = 7 

        thank_you = Text("Thank You!", font_size=60).shift(0.8*UP)
        pietro = Text("Pietro Fumagalli", color= BLUE_C, font_size=35)
        francesco = Text("Francesco Derme", color= BLUE_C, font_size=35)
        pietro.shift(0.15*DOWN).shift(2.2*LEFT)
        francesco.shift(0.15*DOWN).shift(2.2*RIGHT)

        self.wait(1)
        self.play(Write(thank_you))
        self.play(Write(pietro),
                  Write(francesco),
                  run_time = 1.5)
        self.wait(1.5)

        outline_trophies_group = VGroup()
        for _ in range(num_trophies):
            trophy = create_trophy_outline(outline_color=WHITE)
            outline_trophies_group.add(trophy)

        outline_trophies_group.arrange(RIGHT, buff=0.6).scale(0.7)
        outline_trophies_group.next_to(thank_you, 6*DOWN)

        self.play(LaggedStart(
            *[Create(outline) for outline in outline_trophies_group],
            lag_ratio=0.2)
        )

        self.wait(1)

        transform_animations = []
        for i, outline_trophy_instance in enumerate(outline_trophies_group):
            filled_version = create_filled_trophy_shape(fill_color=GOLD)
            filled_version.match_height(outline_trophy_instance)
            filled_version.move_to(outline_trophy_instance.get_center())
            transform_animations.append(Transform(outline_trophy_instance, filled_version))

        self.play(LaggedStart(*transform_animations, lag_ratio=0.2))
        self.wait(2)

        all_2d_elements = VGroup(thank_you, pietro,
                                 francesco, outline_trophies_group)

        self.play(FadeOut(all_2d_elements, shift = 0.5*UP))
        self.wait(0.5)
