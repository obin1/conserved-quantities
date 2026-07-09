from manim import *
import numpy as np


class UnitSimplex3D(ThreeDScene):
    def construct(self):
        # Camera: zoom out a bit so nothing gets cut off at the top
        self.set_camera_orientation(phi=68 * DEGREES, theta=-50 * DEGREES, zoom=0.82)

        # Axes
        axes = ThreeDAxes(
            x_range=[0, 1.35, 0.2],
            y_range=[0, 1.35, 0.2],
            z_range=[0, 1.35, 0.2],
            x_length=6.0,
            y_length=6.0,
            z_length=6.0,
        )

        # Shift everything slightly downward to keep the top inside the frame
        axes.shift(DOWN * 0.45)

        # Axis labels (Text only, no LaTeX needed)
        x_label = Text("RONO2", font_size=28)
        y_label = Text("NO2", font_size=28)
        z_label = Text("NO", font_size=28)

        x_label.move_to(axes.c2p(1.42, 0, 0) + RIGHT * 0.15)
        y_label.move_to(axes.c2p(0, 1.42, 0) + UP * 0.10)
        z_label.move_to(axes.c2p(0, 0, 1.42) + OUT * 0.10)

        # Unit simplex: x + y + z = 1 in the first octant
        p1 = axes.c2p(1, 0, 0)
        p2 = axes.c2p(0, 1, 0)
        p3 = axes.c2p(0, 0, 1)

        simplex = Polygon(
            p1, p2, p3,
            fill_color=GREEN,
            fill_opacity=0.45,
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

        # Second plane from your Plotly code:
        # y = (a1*x + L)/a2
        k1 = 0.8
        k2 = 1.2
        a1 = k1 / (k1 + k2)   # 0.4
        a2 = k2 / (k1 + k2)   # 0.6
        L = 0.1

        def second_plane(u, v):
            y = (a1 * u + L) / a2
            return axes.c2p(u, y, v)

        # A clipped vertical rectangular patch of the plane
        plane2 = Surface(
            second_plane,
            u_range=[0.0, 0.6],
            v_range=[0.0, 0.6],
            resolution=(16, 16),
        )
        plane2.set_fill(PINK, opacity=0.45)
        plane2.set_stroke(PINK, width=0.5)
        plane2.set_shade_in_3d(True)

        # Intersection line from your Plotly code
        t = np.linspace(0, 5 / 6, 300)
        x_line = 0.5 - 3 * t / 5
        y_line = 0.5 - 2 * t / 5
        z_line = t

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

        # Optional title
        title = Text("3D Graph with Unit Simplex", font_size=32)
        title.to_corner(UL)
        title.shift(DOWN * 0.15)

        self.add(axes)
        self.add_fixed_orientation_mobjects(x_label, y_label, z_label, title)

        self.play(FadeIn(simplex))
        self.play(Create(simplex_edges), FadeIn(simplex_vertices))
        self.play(FadeIn(plane2))
        self.play(Create(intersection))

        self.begin_ambient_camera_rotation(rate=0.10)
        self.wait(6)
        self.stop_ambient_camera_rotation()
        self.wait(1)