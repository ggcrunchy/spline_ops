--- Various B&eacute;zier utilities.
--
-- For purposes of this module, an instance of type **Vector** is a value, e.g. a table,
-- that has and / or receives **number** members **x** and **y**.

-- TODO: Investigate
-- "Arc-Length Parameterized Spline Curves for Real-Time Simulation", Hongling Wang, Joseph Kearney, and Kendall Atkinson
-- "Arc Length Parameterization of Spline Curves", John W. Peterson

--
-- Permission is hereby granted, free of charge, to any person obtaining
-- a copy of this software and associated documentation files (the
-- "Software"), to deal in the Software without restriction, including
-- without limitation the rights to use, copy, modify, merge, publish,
-- distribute, sublicense, and/or sell copies of the Software, and to
-- permit persons to whom the Software is furnished to do so, subject to
-- the following conditions:
--
-- The above copyright notice and this permission notice shall be
-- included in all copies or substantial portions of the Software.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
-- EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
-- MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
-- IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
-- CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
-- TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
-- SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
--
-- [ MIT license: http://www.opensource.org/licenses/mit-license.php ]
--

-- Standard library imports --
local abs = math.abs
local log = math.log
local sqrt = math.sqrt

-- Modules --
local arc_length = require("tektite_core.number.arc_length")
local utils = require("spline_ops.utils")

-- Cached module references --
local _GetPosition_
local _GetPositionAndTangent_
local _GetTangent_
local _Length_
local _Length_Array_

-- Exports --
local M = {}

--
--
--

--- Get the position along a quadratic B&eacute;zier spline at time _t_.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q ...control point...
-- @tparam Vector p2 ...and endpoint #2.
-- @number t Interpolation time, &isin; [0, 1].
-- @treturn number Position x-coordinate...
-- @treturn number ...and y-coordinate.
function M.GetPosition (p1, q, p2, t)
	local s = 1 - t
	local a, b, c = s * s, 2 * s * t, t * t

	return a * p1.x + b * q.x + c * p2.x, a * p1.y + b * q.y + c * p2.y
end

--
--
--

--- Array variant of @{GetPosition}.
-- @array bezier Elements 1, 2, 3 are interpreted as arguments _p1_, _q_, _p2_ from @{GetPosition}.
-- @number t Interpolation time, &isin; [0, 1].
-- @treturn number Position x-coordinate...
-- @treturn number ...and y-coordinate.
function M.GetPosition_Array (bezier, t)
	return _GetPosition_(bezier[1], bezier[2], bezier[3], t)
end

--
--
--

--- Variant of @{GetPosition} that also returns the tangent.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q ...control point...
-- @tparam Vector p2 ...and endpoint #2.
-- @number t Interpolation time, &isin; [0, 1].
-- @treturn number Position x-coordinate...
-- @treturn number ...and y-coordinate.
-- @treturn number Tangent x-component...
-- @treturn number ...and y-component.
function M.GetPositionAndTangent (p1, q, p2, t)
	local s, t2 = 1 - t, 2 * t
	local a, b, c = s * s, 2 * s * t, t * t
	local px, py = a * p1.x + b * q.x + c * p2.x, a * p1.y + b * q.y + c * p2.y
	local d, e, f = t2 - 2, 2 * (1 - t2), t2
	local tx, ty = d * p1.x + e * q.x + f * p2.x, d * p1.y + e * q.y + f * p2.y

	return px, py, tx, ty
end

--
--
--

--- Array variant of @{GetPositionAndTangent}.
-- @array bezier Elements 1, 2, 3 are interpreted as arguments _p1_, _q_, _p2_ from @{GetPosition}.
-- @number t Interpolation time, &isin; [0, 1].
-- @treturn number Position x-coordinate...
-- @treturn number ...and y-coordinate.
-- @treturn number Tangent x-component...
-- @treturn number ...and y-component.
function M.GetPositionAndTangent_Array (bezier, t)
	return _GetPositionAndTangent_(bezier[1], bezier[2], bezier[3], t)
end

--
--
--

--- Variant of @{GetPosition} that returns the tangent.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q ...control point...
-- @tparam Vector p2 ...and endpoint #2.
-- @number t Interpolation time, &isin; [0, 1].
-- @treturn number Tangent x-component...
-- @treturn number ...and y-component.
function M.GetTangent (p1, q, p2, t)
	local t2 = 2 * t
	local a, b, c = t2 - 2, 2 * (1 - t2), t2

	return a * p1.x + b * q.x + c * p2.x, a * p1.y + b * q.y + c * p2.y
end

--
--
--

--- Array variant of @{GetTangent}.
-- @array bezier Elements 1, 2, 3 are interpreted as arguments _p1_, _q_, _p2_ from @{GetPosition}.
-- @number t Interpolation time, &isin; [0, 1].
-- @treturn number Tangent x-component...
-- @treturn number ...and y-component.
function M.GetTangent_Array (bezier, t)
	return _GetTangent_(bezier[1], bezier[2], bezier[3], t)
end

--
--
--

-- Gets control point (as perpendicular displacement from midpoint)
local function AuxGetControlPoint (x1, y1, x2, y2, below)
	local midx, midy = x1 + .5 * (x2 - x1), y1 + .5 * (y2 - y1)
	local to_x, to_y = midx - x1, midy - y1

	if below then
		return midx - to_y, midy + to_x
	else
		return midx + to_y, midy - to_x
	end
end

--- Compute a reasonable control point for a quadratic B&eacute;zier spline.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector p2 ...and #2.
-- @bool below Should the control points be "below" the line segment between _p1_ and
-- _p2_? (For this purpose, _p1_ is considered to be on the left, _p2_ on the right.)
-- @treturn number Position x-coordinate...
-- @treturn number ...and y-coordinate.
function M.GetControlPoint (p1, p2, below)
	return AuxGetControlPoint(p1.x, p1.y, p2.x, p2.y, below)
end

--
--
--

--- Compute a reasonable control point for a quadratic B&eacute;zier spline.
--
-- When the endpoints do not line up (horizontally or vertically), they may be interpreted
-- as two opposite corners of a rectangle, and one of the unused corners is chosen as the
-- control point.
--
-- Otherwise, the algorithm can either fall back to the behavior of @{GetControlPoint} or
-- choose the midpoint of _p1_ and _p2_ (i.e. the spline degenerates to a line).
--
-- For purposes of above / below, cf. _below_ in @{GetControlPoint}.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector p2 ...and #2.
-- @number too_close If the difference (rather, its absolute value) between the x- or
-- y-coordinates of _p1_ and _p2_ is less than this amount, the points are considered
-- as lined up, and the corner is abandoned.
-- @string[opt] how If a corner can be found, the control point "below" the segment is chosen
-- if _how_ is **"below"** or **"below\_else\_middle"**, or the one "above" otherwise.
--
-- Failing that, the midpoint is chosen as a fallback if _how_ is **"above\_else\_middle"**
-- or **below\_else\_middle"**. Otherwise, @{GetControlPoint} is the fallback, with _below_
-- true if _how_ is **"below"**.
-- @treturn number Position x-coordinate...
-- @treturn number ...and y-coordinate.
function M.GetControlPoint_TryCorner (p1, p2, too_close, how)
	local x1, y1, x2, y2 = p1.x, p1.y, p2.x, p2.y

	if abs(x2 - x1) > too_close and abs(y2 - y1) > too_close then
		-- TODO: Hack, reason the above/below out properly
		if x2 < x1 then
			x1, y1, x2, y2 = x2, y2, x1, y1
		end

		-- Choose one of the corners.
		local ax, ay, bx, by = x1, y2, x2, y1

		if y2 < y1 then
			ax, ay, bx, by = bx, by, ax, ay
		end

		if how == "below" or "below_else_middle" then
			return bx, by
		else
			return ax, ay
		end
	end

	-- Degenerate rectangle (no free corners): fall back to something else.
	if how == "above_else_middle" or how == "below_else_middle" then
		return .5 * (x1 + x2), .5 * (y1 + y2)
	else
		return AuxGetControlPoint(x1, y1, x2, y2, how == "below")
	end
end

--
--
--

--- Compute a (degree 2) [B&eacute;zier spline's length](http://malczak.info/blog/quadratic-bezier-curve-length/).
-- @tparam Vector p1 Endpoint #1 of control polygon...
-- @tparam Vector q ...interior control point...
-- @tparam Vector p2 ...and endpoint #2.
-- @treturn number Approximate arc length.
function M.Length (p1, q, p2)
	local p1x, p1y = p1.x, p1.y
	local qpx, qpy = 2 * q.x - p1x, 2 * q.y - p1y
	local ax, bx = p2.x - qpx, qpx - p1x
	local ay, by = p2.y - qpy, qpy - p1y

	local A = ax * ax + ay * ay
	local C = bx * bx + by * by

	if A > 1e-9 then
		A = 4 * A

		local B = 4 * (ax * bx + ay * by)
		local Sabc = 2 * sqrt(A + B + C)
		local A_2, C_2 = sqrt(A), 2 * sqrt(C)
		local A_32, BA = 2 * A * A_2, B / A_2

		return (A_32 * Sabc + A_2 * B * (Sabc - C_2) + (4 * C * A - B * B) * log((2 * A_2 + BA + Sabc) / (BA + C_2))) / (4 * A_32)
	else
		return sqrt(C)
	end
end

--
--
--

--- Array variant of @{Length}.
-- @array bezier Elements 1, 2, 3 are interpreted as arguments _p1_, _q_, _p2_ from @{Length}.
-- @treturn number Approximate arc length.
function M.Length_Array (bezier)
	return _Length_(bezier[1], bezier[2], bezier[3])
end

--
--
--

local function Subdivide (bezier, left, right, t)
	return utils.Subdivide(bezier, left, right, t, 2)
end

local Bezier = {}

--- Populate an arc &rarr; parameter lookup table given a (degree 2) B&eacute;zier spline.
-- @array lut Lookup table, cf. @{tektite_core.number.sampling.Lookup}.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q ...control point...
-- @tparam Vector p2 ...and endpoint #2.
-- @int[opt] nsamples Number of samples to load into _lut_. If absent, a default is used.
-- @treturn number Total arc length.
function M.PopulateArcLengthLUT (lut, p1, q, p2, nsamples)
	Bezier[1], Bezier[2], Bezier[3] = p1, q, p2

	local s = arc_length.PopulateLUT_Subdivide(lut, Bezier, Subdivide, nsamples, _Length_Array_)

	Bezier[1], Bezier[2], Bezier[3] = nil

	return s
end

--
--
--

_GetPosition_ = M.GetPosition
_GetPositionAndTangent_ = M.GetPositionAndTangent
_GetTangent_ = M.GetTangent
_Length_ = M.Length
_Length_Array_ = M.Length_Array

return M