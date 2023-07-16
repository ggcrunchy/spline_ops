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
local ipairs = ipairs

-- Modules --
local arc_length = require("tektite_core.number.arc_length")
local utils = require("spline_ops.utils")

-- Extension imports --
local hypot = math.hypot

-- Cached module references --
local _Length_Array_

-- Exports --
local M = {}

--
--
--

--[[
	Length via split method

	Earl Boebert's commentary, in original:

	The last suggestion by Gravesen is pretty nifty, and I think it's a candidate for the
	next Graphics Gems. I hacked out the following quick implementation, using the .h and
	libraries definitions from Graphics Gems I (If you haven't got that book then you have
	no business mucking with with this stuff :-)) The function "bezsplit" is lifted
	shamelessly from Schneider's Bezier curve-fitter.
]]

-- Workspace matrix: a triangle of vectors is used; the rest is dummied out --
local Temp = {}

for i = 1, 4 do
	Temp[i] = {}

	for j = 1, 4 do
		Temp[i][j] = j <= 5 - i and {}
	end
end

-- Push a point onto the stack
local V, Top = {}

local function AddPoint (point)
	V[Top + 1], V[Top + 2], Top = point.x, point.y, Top + 2
end

-- Split a cubic bezier in two
local function BezSplit ()
	-- Copy control points.
	local base = Top + 1

	for _, temp in ipairs(Temp[1]) do
		temp.x, temp.y, base = V[base], V[base + 1], base + 2
	end

	-- Triangle computation.
	local prev_row = Temp[1]

	for i = 2, 4 do
		local row = Temp[i]

		for j = 1, 5 - i do
			local r, pr1, pr2 = row[j], prev_row[j], prev_row[j + 1]

			r.x, r.y = .5 * (pr1.x + pr2.x), .5 * (pr1.y + pr2.y)
		end

		prev_row = row
	end

	-- Left split.
	for i = 1, 4 do
		AddPoint(Temp[i][1])
	end

	-- Right split.
	for i = 1, 4 do
		AddPoint(Temp[5 - i][i])
	end
end

-- Add polyline length if close enough
local function AddIfClose (length, err)
	-- Pop four points off the stack. Compute the point-to-point and chord lengths.
	Top = Top - 8

	local base = Top + 1
	local x, y = V[base], V[base + 1]

	local len, main_len = 0, hypot(V[base + 6] - x, V[base + 7] - y)

	for _ = 1, 3 do
		local dx, dy = V[base + 2] - x, V[base + 3] - y

		len = len + hypot(dx, dy)
		base, x, y = base + 2, x + dx, y + dy
	end

	-- If the point-to-point lengths sum to much more than the chord length, split up
	-- the curve and sum the lengths of the two parts.
	if len - main_len > err then
		BezSplit()

		local ll = AddIfClose(length, err)
		local lr = AddIfClose(length, err)

		len = ll + lr
	end

	return len
end

--- Compute a (degree 3) [B&eacute;zier spline's length](http://steve.hollasch.net/cgindex/curves/cbezarclen.html).
-- @tparam Vector p1 Endpoint #1 of control polygon...
-- @tparam Vector q1 ...interior control point #1...
-- @tparam Vector q2 ...interior control point #2...
-- @tparam Vector p2 ...and endpoint #2.
-- @number tolerance "Close enough" separation between arc estimates and chord lengths.
-- @treturn number Approximate arc length.
function M.Length (p1, q1, q2, p2, tolerance)
	Top = 0

	AddPoint(p1)
	AddPoint(q1)
	AddPoint(q2)
	AddPoint(p2)

	return AddIfClose(0, tolerance)
end

--
--
--

--- Array variant of @{Length3}.
-- @array bezier Elements 1, 2, 3, 4 are interpreted as arguments _p1_, _q1_, _q2_, _p2_
-- from @{Length3}.
-- @number tolerance "Close enough" separation between arc estimates and chord lengths.
-- @treturn number Approximate arc length.
function M.Length_Array (bezier, tolerance)
	Top = 0

	for i = 1, 4 do
		AddPoint(bezier[i])
	end

	return AddIfClose(0, tolerance)
end

--
--
--

local function Subdivide (bezier, left, right, t)
	return utils.Subdivide(bezier, left, right, t, 3)
end

local Bezier = {}

--- Populate an arc &rarr; parameter lookup table given a (degree 3) B&eacute;zier spline.
-- @array lut Lookup table, cf. @{tektite_core.number.sampling.Lookup}.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q1 ...control point #1...
-- @tparam Vector q2 ...control point #2...
-- @tparam Vector p2 ...and endpoint #2.
-- @int[opt] nsamples Number of samples to load into _lut_. If absent, a default is used.
-- @number[opt] tolerance "Close enough" tolerance, cf. @{spline_ops.bezier3.Length}.
-- If absent, a default is used.
-- @treturn number Total arc length.
function M.PopulateArcLengthLUT (lut, p1, q1, q2, p2, nsamples, tolerance)
	Bezier[1], Bezier[2], Bezier[3], Bezier[4] = p1, q1, q2, p2

	local s = arc_length.PopulateLUT_Subdivide(lut, Bezier, Subdivide, nsamples, _Length_Array_, tolerance or 1e-3)

	Bezier[1], Bezier[2], Bezier[3], Bezier[4] = nil

	return s
end

--
--
--

_Length_Array_ = M.Length_Array

return M