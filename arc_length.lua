--- Some utilities for building arc-length lookup tables, built atop @{tektite_core.number.sampling}.
--
-- The samples' **x** values correspond to the arc length, _s_, and go from 0 to the full
-- length of the curve. The **y** values correspond to a "curve time", _t_, that increases
-- monotonically from 0 to 1.
--
-- For purposes of this module, an instance of type **Vector** is a value, e.g. a table,
-- that has and / or receives **number** members **x** and **y**. (**N.B.** The samples from
-- the previous paragraph are unrelated, despite the coincidence of field names.)

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

-- Modules --
local bezier = require("spline_ops.bezier")
local integrators = require("tektite_core.number.integrators")
local sampling = require("tektite_core.number.sampling")

-- Exports --
local M = {}

-- Adds the final "full" arc length to the LUT and makes it ready to use
local function CloseLUT (lut, s)
	sampling.AddSample(lut, s, 1)

	return s
end

-- Intermediate spline vectors, partitions --
local Bezier, Left, Right = {}, {}, {}

-- Build up the LUT from a (degree n) B&eacute;zier spline
local function SetLUT_Bezier (lut, nsamples, func, tolerance)
	nsamples = nsamples or 20

	local spline, deg = Bezier, #Bezier - 1
	local s, t, dt = 0, 0, 1 / nsamples

	repeat
		sampling.AddSample(lut, s, t)

		-- Divide the curve into parts of length u = 1 / nsamples. On the first iteration,
		-- the subdivision parameter is trivially u itself, leaving a right-hand side of
		-- length (nsamples - 1) / nsamples. On the second iteration, to maintain length u,
		-- we have 1 / nsamples = t * (nsamples - 1) / nsamples, i.e. new parameter t = 1 /
		-- (nsamples - 1). In general, on interation i, t = 1 / (nsamples - i + 1). (On the
		-- final iteration, t = 1, and the right-hand side is empty.)
		bezier.Subdivide(spline, Left, Right, 1 / nsamples, deg)

		local ds = func(Left, tolerance)

		spline, s, t, nsamples = Right, s + ds, t + dt, nsamples - 1
	until nsamples == 0

	return CloseLUT(lut, s)
end

--- Populates an arc &rarr; parameter lookup table given a (degree 2) B&eacute;zier spline.
-- @array lut Lookup table, cf. @{Lookup}.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q ...control point...
-- @tparam Vector p2 ...and endpoint #2.
-- @int[opt] nsamples Number of samples to load into _lut_. If absent, a default is used.
-- @treturn number Total arc length.
function M.SetLUT_Bezier2 (lut, p1, q, p2, nsamples)
	Bezier[1], Bezier[2], Bezier[3] = p1, q, p2

	local s = SetLUT_Bezier(lut, nsamples, bezier.Length2_Array)

	Bezier[1], Bezier[2], Bezier[3] = nil

	return s
end

--- Populates an arc &rarr; parameter lookup table given a (degree 3) B&eacute;zier spline.
-- @array lut Lookup table, cf. @{Lookup}.
-- @tparam Vector p1 Endpoint #1...
-- @tparam Vector q1 ...control point #1...
-- @tparam Vector q2 ...control point #2...
-- @tparam Vector p2 ...and endpoint #2.
-- @int[opt] nsamples Number of samples to load into _lut_. If absent, a default is used.
-- @number[opt] tolerance "Close enough" tolerance, cf. @{spline_ops.bezier.Length3}.
-- If absent, a default is used.
-- @treturn number Total arc length.
function M.SetLUT_Bezier3 (lut, p1, q1, q2, p2, nsamples, tolerance)
	Bezier[1], Bezier[2], Bezier[3], Bezier[4] = p1, q1, q2, p2

	local s = SetLUT_Bezier(lut, nsamples, bezier.Length3_Array, tolerance or 1e-3)

	Bezier[1], Bezier[2], Bezier[3], Bezier[4] = nil

	return s
end

--- Populates an arc &rarr; parameter lookup table given a function to integrate over [0, 1].
-- @array lut Lookup table, cf. @{Lookup}.
-- @string? how If this is **"gauss_legendre"**, @{tektite_core.number.integrators.GaussLegendre} is used
-- as the integration method. Otherwise, @{tektite_core.number.integrators.Romberg} is used.
-- @callable func Function to integrate, e.g. an integrand supplied by @{spline_ops.cubic.LineIntegrand}.
-- @int[opt] nsamples Number of samples to load into _lut_. If absent, a default is used.
-- @number[opt] tolerance Tolerance, as used by some integrators. If absent, a default is used.
-- @treturn number Total arc length.
function M.SetLUT_Func (lut, how, func, nsamples, tolerance)
	nsamples, tolerance = nsamples or 20, tolerance or 1e-3

	if how == "gauss_legendre" then
		how = integrators.GaussLegendre
	else
		how = integrators.Romberg
	end

	local a, s, dt = 0, 0, 1 / nsamples

	for _ = 1, nsamples do
		sampling.AddSample(lut, s, a)

		local b = a + dt
		local ds = how(func, a, b, tolerance)

		a, s = b, s + ds
	end

	return CloseLUT(lut, s)
end

-- Export the module.
return M