module SimpleGeometry

export Point, Line, Circle, vector, intersection, rotate

using StaticArrays

struct Point
    xy::SVector{2, Float64}
end

magnitude(p::Point) = sqrt(sum(p.xy .^ 2))
normalize(p::Point) = p / magnitude(p)

Point(x, y) = Point(SVector{2, Float64}(x, y))

function Base.getproperty(p::Point, sym::Symbol)
    if sym == :x
        return p.xy[1]
    elseif sym == :y
        return p.xy[2]
    else
        getfield(p, sym)
    end
end

Base.:+(p1::Point, p2::Point) = Point(p1.xy + p2.xy)
Base.:-(p1::Point, p2::Point) = Point(p1.xy - p2.xy)
Base.:*(p::Point, factor::Real) = Point(p.xy .* factor)
Base.:*(factor::Real, p::Point) = p * factor
Base.:/(p::Point, r::Real) = Point(p.xy ./ r)
from_to(p1::Point, p2::Point) = p2 - p1
between(p1::Point, p2::Point, fraction::Real) = p1 + (p2 - p1) * fraction
cross(p1::Point, p2::Point) = p1.x * p2.y - p1.y * p2.x
dot(p1::Point, p2::Point) = p1.x * p2.x + p1.y * p2.y

function angle(p::Point; degrees=True)
    radians = atan(p.y, p.x)
    degrees ? rad2deg(radians) : radians
end

function signed_angle_to(p1::Point, p2::Point; degrees::Bool=true)
    radians = atan(cross(p1, p2), dot(p1, p2))
    return degrees ? rad2deg(radians) : radians
end

function rotate(p::Point, angle::Real; around::Point=Point(0, 0), degrees=true)
    vector = from_to(around, p)
    rotated_vector = Point(rotation_matrix(angle, degrees=degrees) * vector.xy)
    rotated_vector + around
end

function rotation_matrix(angle; degrees=true)
    angle = degrees ? deg2rad(angle) : angle
    c = cos(angle)
    s = sin(angle)
    SMatrix{2, 2, Float64}(c, s, -s, c)
end

struct Line
    from::Point
    to::Point
end

xs(l::Line) = SVector(l.from.x, l.to.x)
ys(l::Line) = SVector(l.from.y, l.to.y)

function intersection(l1::Line, l2::Line)
    x1, x2 = xs(l1)
    y1, y2 = ys(l1)
    x3, x4 = xs(l2)
    y3, y4 = ys(l2)
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) /
        ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))

    px = x1 + t * (x2 - x1)
    py = y1 + t * (y2 - y1)
    Point(px, py)
end

vector(l::Line) = from_to(l.from, l.to)
angle(l::Line) = angle(vector(l))
length(l::Line) = magnitude(vector(l))
fraction(l::Line, frac::Real) = between(l.from, l.to, frac)
reversed(l::Line) = Line(l.to, l.from)
direction(l::Line) = normalize(vector(l))

function scale(l::Line, scalar::Real)
    movement = (scalar - 1) * vector(l) / 2
    Line(l.from - movement, l.to + movement)
end

function scaleto(l::Line, len::Real)
    scalar = len / length(l)
    scale(l, scalar)
end

function addlength(l::Line, len::Real)
    dir = direction(l)
    movement = len / 2 * dir
    Line(l.from - movement, l.to + movement)
end

function rotate(l::Line, angle::Real; around::Point=Point(0, 0), degrees=true)
    Line(
        rotate(l.from, angle, around=around, degrees=degrees),
        rotate(l.to, angle, around=around, degrees=degrees)
    )
end

struct Circle
    center::Point
    radius::Float64
end

function Circle(p1::Point, p2::Point, p3::Point)
    l1 = Line(p1, p2)
    c1 = fraction(l1, 0.5)
    l2 = Line(p2, p3)
    c2 = fraction(l2, 0.5)

    dl1 = vector(l1)
    dl2 = vector(l2)
    cl1 = Line(c1, c1 + Point(dl1.y, -dl1.x))
    cl2 = Line(c2, c2 + Point(dl2.y, -dl2.x))

    center = intersection(cl1, cl2)
    radius = from_to(p1, center) |> magnitude
    Circle(center, radius)
end

function Circle(center::Point, p1::Point)
    radius = from_to(center, p1) |> magnitude
    Circle(center, radius)
end

function intersection(c::Circle, l::Line)
    # algorithm works for circle at (0, 0)
    frm_moved = l.from - c.center
    to_moved = l.to - c.center
    dx, dy = from_to(frm_moved, to_moved).xy
    dr = sqrt(dx ^ 2 + dy ^ 2)
    D = cross(frm_moved, to_moved)

    delta = (c.radius ^ 2) * (dr ^ 2) - (D ^ 2)

    if abs(delta) < 1e-10
        # tangent
        x = D * dy / (dr ^ 2)
        y = -D * dx / (dr ^ 2)
        return Point(x, y) + c.center
    elseif delta < 0
        return nothing
    else
        xplusminus = sign(dy) * dx * sqrt(delta)
        yplusminus = abs(dy) * sqrt(delta)

        x1 = (D * dy + xplusminus) / (dr ^ 2)
        x2 = (D * dy - xplusminus) / (dr ^ 2)
        y1 = (-D * dx + yplusminus) / (dr ^ 2)
        y2 = (-D * dx - yplusminus) / (dr ^ 2)
        return (Point(x1, y1) + c.center, Point(x2, y2) + c.center)
    end
end

function tangentpoints(c::Circle, through::Point)
    p_to_center = from_to(through, c.center)
    a = asin(c.radius / magnitude(p_to_center))
    b = atan(p_to_center.y, p_to_center.x)
    t1 = b - a
    p_tangent_1 = Point(sin(t1), -cos(t1)) * c.radius + c.center
    t2 = b + a
    p_tangent_2 = Point(-sin(t2), cos(t2)) * c.radius + c.center
    return p_tangent_1, p_tangent_2
end

function tangents(c::Circle, through::Point)
    p_tangent_1, p_tangent_2 = tangentpoints(c, through)
    return Line(through, p_tangent_1), Line(through, p_tangent_2)
end

function outertangents(c1::Circle, c2::Circle)
    big_circle, small_circle = (c1, c2) if c1.radius > c2.radius else (c2, c1)
    radius_difference = big_circle.radius - small_circle.radius
    gamma = -atan(big_circle.center.y - small_circle.center.y, big_circle.center.x - small_circle.center.x)

    beta = asin(radius_difference / magnitude(from_to(small_circle.center, big_circle.center)))
    alpha1 = gamma - beta
    alpha2 = gamma + beta

    small_tangent_point_1 = small_circle.center + Point(
        small_circle.radius * cos(pi / 2 - alpha1),
        small_circle.radius * sin(pi / 2 - alpha1),
    )

    big_tangent_point_1 = big_circle.center + Point(
        big_circle.radius * cos(pi / 2 - alpha1),
        big_circle.radius * sin(pi / 2 - alpha1),
    )

    small_tangent_point_2 = small_circle.center + Point(
        small_circle.radius * cos(-pi / 2 - alpha2),
        small_circle.radius * sin(-pi / 2 - alpha2),
    )

    big_tangent_point_2 = big_circle.center + Point(
        big_circle.radius * cos(-pi / 2 - alpha2),
        big_circle.radius * sin(-pi / 2 - alpha2),
    )

    return Line(small_tangent_point_1, big_tangent_point_1), Line(small_tangent_point_2, big_tangent_point_2)
end

function point_at_angle(c::Circle, angle::Real; degrees::Bool=true)
    angle = degrees ? deg2rad(degrees) : angle
    c.center + Point(cos(angle), sin(angle)) * c.radius
end

function closestto(c::Circle, p::Point)
    normalize(from_to(c.center, p)) * c.radius
end

area(c::Circle) = pi * (c.radius ^ 2)
circumference(c::Circle) = 2pi * r

struct CircularArc
    center::Point
    radius::Float64
    start_angle::Float64
    end_angle::Float64
    degrees::Bool
end

function CircularArc(p1::Point, p2::Point, p3::Point, degrees::Bool)
    circle = Circle(p1, p2, p3)
    start_angle = angle(from_to(circle.center, p1), degrees=degrees)
    end_angle = angle(from_to(circle.center, p3), degrees=degrees)

    if signed_angle_to(from_to(p1, p3), from_to(p1, p2)) > 0
        start_angle, end_angle = end_angle, start_angle
    end
    CircularArc(circle.center, circle.radius, start_angle, end_angle, degrees)
end

struct Rect
    center::Point
    width::Float64
    height::Float64
    angle::Float64
    degrees::Bool
end

function lowerleft(r::Rect)
    r.center + rotate(Point(-r.width * 0.5, -r.height * 0.5), r.angle, degrees=r.degrees)
end

function lowerright(r::Rect)
    r.center + rotate(Point( r.width * 0.5, -r.height * 0.5), r.angle, degrees=r.degrees)
end

function upperleft(r::Rect)
    r.center + rotate(Point(-r.width * 0.5,  r.height * 0.5), r.angle, degrees=r.degrees)
end

function upperright(r::Rect)
    r.center + rotate(Point( r.width * 0.5,  r.height * 0.5), r.angle, degrees=r.degrees)
end

topline(r::Rect) = Line(upperleft(r), upperright(r))
bottomline(r::Rect) = Line(lowerleft(r), lowerright(r))
leftline(r::Rect) = Line(lowerleft(r), upperleft(r))
rightline(r::Rect) = Line(lowerright(r), upperright(r))

struct RoundRect
end

end # module
