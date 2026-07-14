from manim import *
import numpy as np

class UnitSimplex3D(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=68 * DEGREES, theta=-50 * DEGREES, zoom=1.05)

        axes = ThreeDAxes(
            x_range=[0, 1.35, 0.2],
            y_range=[0, 1.35, 0.2],
            z_range=[0, 1.35, 0.2],
            x_length=3.6,
            y_length=3.6,
            z_length=3.6,
            x_axis_config={"include_ticks": False, "include_numbers": False},
            y_axis_config={"include_ticks": False, "include_numbers": False},
            z_axis_config={"include_ticks": False, "include_numbers": False},
        )
        axes.shift(DOWN * 2.5)
        axes.shift(RIGHT * 1.5)

        x_label = Text("RONO₂", font_size=28)
        y_label = Text("NO₂", font_size=28)
        z_label = Text("NO", font_size=28)

        x_label.move_to(axes.c2p(1.42, 0, 0)).shift(RIGHT * 0.25)
        y_label.move_to(axes.c2p(0, 1.42, 0)).shift(OUT * 0.25)
        z_label.move_to(axes.c2p(0, 0, 1.42)).shift(RIGHT * 0.10)

        # Unit simplex x + y + z = 1
        p1 = axes.c2p(1, 0, 0)
        p2 = axes.c2p(0, 1, 0)
        p3 = axes.c2p(0, 0, 1)

        simplex = Polygon(
            p1, p2, p3,
            fill_color=GREEN,
            fill_opacity=0.75,
            stroke_color=GREEN,
            stroke_width=2,
        )
        simplex.set_shade_in_3d(True)

        simplex_edges = VGroup(
            Line3D(p1, p2, color=GREEN, thickness=0.03),
            Line3D(p2, p3, color=GREEN, thickness=0.03),
            Line3D(p3, p1, color=GREEN, thickness=0.03),
        )

        simplex_vertices = VGroup(
            Dot3D(p1, color=WHITE, radius=0.05),
            Dot3D(p2, color=WHITE, radius=0.05),
            Dot3D(p3, color=WHITE, radius=0.05),
        )

        # Contour-like lines on the simplex: x = const, y = const, z = const
        simplex_contours = VGroup()

        # x = c lines on x + y + z = 1
        for c in np.linspace(0, 1, 20):
            simplex_contours.add(
                Line3D(
                    axes.c2p(c, 0, 1 - c),
                    axes.c2p(c, 1 - c, 0),
                    color=GREEN,
                    thickness=0.010,
                )
            )

        # y = c lines
        for c in np.linspace(0, 1, 20):
            simplex_contours.add(
                Line3D(
                    axes.c2p(0, c, 1 - c),
                    axes.c2p(1 - c, c, 0),
                    color=GREEN,
                    thickness=0.010,
                )
            )

        # z = c lines
        for c in np.linspace(0, 1, 20):
            simplex_contours.add(
                Line3D(
                    axes.c2p(0, 1 - c, c),
                    axes.c2p(1 - c, 0, c),
                    color=GREEN,
                    thickness=0.010,
                )
            )

        # y = (a1*x + L)/a2
        k1 = 0.8
        k2 = 1.2
        a1 = k1 / (k1 + k2)   # 0.4
        a2 = k2 / (k1 + k2)   # 0.6
        L = 0.1

        # def second_plane(u, v):
        #     y = (a1 * u + L) / a2
        #     return axes.c2p(u, y, v)

        # plane2 = Surface(
        #     second_plane,
        #     u_range=[0.0, 0.60],   # further in the y direction
        #     v_range=[0.0, 1.00],   # higher in z
        #     resolution=(24, 24),
        # )
        # plane2.set_fill(PINK, opacity=1.0)
        # plane2.set_stroke(PINK, width=0.5)
        # plane2.set_shade_in_3d(True)

        def plane_left(u, t):
            v_line = 5/6 - 5*u/3
            v = t * v_line            # 0 <= v <= intersection
            y = (a1*u + L)/a2
            return axes.c2p(u, y, v)

        plane_left = Surface(
            plane_left,
            u_range=[0, 0.5],
            v_range=[0, 1],
        )

        plane_left.set_fill(PINK, opacity=1.0)
        plane_left.set_stroke(PINK, width=0.5)
        plane_left.set_shade_in_3d(True)

        def plane_right(u, t):
            v_line = 5/6 - 5*u/3
            v = v_line + t*(1 - v_line)
            y = (a1*u + L)/a2
            return axes.c2p(u, y, v)

        plane_right = Surface(
            plane_right,
            u_range=[0, 0.5],
            v_range=[0, 1],
        )

        plane_right.set_fill(PINK, opacity=1.0)
        plane_right.set_stroke(PINK, width=0.5)
        plane_right.set_shade_in_3d(True)  

        # Intersection line
        intersection = ParametricFunction(
            lambda s: axes.c2p(
                0.5 - 3 * s / 5,
                0.5 - 2 * s / 5,
                s
            ),
            t_range=[0, 5 / 6],
            color=WHITE,
            stroke_width=6,
        )
        

        # Arrow annotation to the intersection line
        # Arrow annotation to the intersection line
        # arrow = Arrow3D(
        #     start=axes.c2p(1.05, 1.05, 1.05),
        #     end=axes.c2p(0.34, 0.39, 0.45),
        #     color=WHITE,
        # )

        # IMPORTANT: disable depth test on every submobject of the arrow
        # for mob in arrow.family_members_with_points():
        #     mob.set_depth_test(False)
        #     mob.set_z_index(1000)

        # Make the arrow a bit thicker so it reads clearly
        # arrow.set_stroke(width=8)

        # arrow.set_depth_test(False)
        # arrow.set_z_index(100)
        
        # line_label = Text("intersection line", font_size=24)
        # line_label.move_to(axes.c2p(1.75, 1.45, 1.20)).shift(RIGHT * 0.6 + UP * 0.3)
        # line_label.set_color(WHITE)
        # line_label.set_depth_test(False)

        # self.add_fixed_orientation_mobjects(line_label)

        # Add everything
        self.add(axes)
        self.add_fixed_orientation_mobjects(x_label, y_label, z_label)
        self.play(FadeIn(simplex))
        self.play(Create(simplex_edges), FadeIn(simplex_vertices))
        self.play(Create(simplex_contours))
        # self.play(FadeIn(plane2))
        self.play(FadeIn(plane_left))
        self.play(FadeIn(plane_right))
        self.play(Create(intersection))
        # self.play(Create(arrow), FadeIn(line_label))
        # self.add(arrow)
        # self.add(line_label)

        self.begin_ambient_camera_rotation(rate=0.10)
        self.wait(6)
        self.stop_ambient_camera_rotation()
        self.wait(1)