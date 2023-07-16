--- Adapted from **quangnle**'s implementation of P.J. Schneider's [curve-fitting algorithm](https://www.realtimerendering.com/resources/GraphicsGems/gems/FitCurves.c),
-- as found in [CurveFitting](https://github.com/quangnle/CurveFitting_PJSchneider).
--
-- For purposes of this module, an instance of type **Vector** is a value, e.g. a table,
-- that has and / or receives **number** members **x** and **y**.

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

-- Extension imports --
local hypot = math.hypot
local normalize = math.normalize

-- Exports --
local M = {}

--
--
--

local function Bernstein (t)
	local s = 1 - t
	local _3st = 3 * s * t

	return s^3, _3st * s, _3st * t, t^3
end

--
--
--

local Temp = {}

--
--
--

local function ComputeParam (degree, points, t)
	-- Copy control points.
	local n = 0

	for i = 1, degree + 1 do
		local p = points[i]

		n, Temp[n + 1], Temp[n + 2] = n + 2, p.x, p.y		
	end
	
	-- Triangle computation.
	for i = 2, degree + 1 do
		local index, x1, y1 = 0, Temp[1], Temp[2]

		for _ = 1, degree - i do
			local x2, y2 = Temp[index + 3], Temp[index + 4]

			Temp[index], Temp[index + 1] = x1 + t * (x2 - x1), y1 + t * (y2 - y1)

			index, x1, y1 = index + 2, x2, y2
		end
	end

	return Temp[1], Temp[2]
end

--
--
--

local Q1, Q2 = {}, {}

--
--
--

local function GenerateControlPoints (to, from, n)
	local p1 = from[1]
	local p1x, p1y = p1.x, p1.y

	for i = 1, n do
		local p2, out = from[i + 1], to[i]
		local p2x, p2y = p2.x, p2.y
		
		out.x = (p2x - p1x) * n -- n.b. scale factor and loop count happen to coincide in each case
		out.y = (p2y - p1y) * n
		p1x, p1y = p2x, p2y
	end
end

--
--
--

local function FindRoot (points, p, u)
	-- Compute Q'(u) and Q''(u).
	local qux, quy = ComputeParam(3, points, u)

	GenerateControlPoints(Q1, points, 3)
	GenerateControlPoints(Q2, Q1, 2)

	local qu1x, qu1y = ComputeParam(2, Q1, u)
	local qu2x, qu2y = ComputeParam(1, Q2, u)
	
	-- Compute f(u) / f'(u).
	local dx, dy = qux - p.x, quy - p.y
	local numer = dx * qu1x + dy * qu1y
	local denom = qu1x * qu1x + qu1y * qu1y + dx * qu2x + dy * qu2y
	
	if 1 + denom * denom == 1 then -- root found already
		return u
	end

	-- Newton: t <- t - f(x) / f'(x)
	return u - numer / denom
end

--
--
--

-- to solve to matrix equation
-- | C11 C12 | |A1| = |X1|
-- | C21 C22 | |A2|   |X2|
-- A1, A2 are alpha1, alpha2
local function Generate (points, i1, i2, uPrime, tx1, ty1, tx2, ty2)
	local vf, vl, j = points[i1], points[i2], 1
	local vfx, vfy = vf.x, vf.y
	local vlx, vly = vl.x, vl.y

	--
	--
	--

	local C11, C12, C22, X1, X2 = 0, 0, 0, 0, 0
	
	for i = i1, i2 do
		local a, b, c, d = Bernstein(uPrime[j])

		-- compute A's (alpha1 and alpha2)
		local A1x, A1y = tx1 * b, ty1 * b
		local A2x, A2y = tx2 * c, ty2 * c
		
		-- update X
		C11 = C11 + A1x * A1x + A1y * A1y
		C12 = C12 + A1x * A2x + A1y * A2y
		C22 = C22 + A2x * A2x + A2y * A2y
		
		local vfi = points[i]
		local x12, y12 = vfx * (a + b), vfy * (a + b)
		local x34, y34 = vlx * (c + d), vly * (c + d)
		local sumx, sumy = vfi.x - (x12 + x34), vfi.y - (y12 + y34)
		
		X1 = X1 + A1x * sumx + A1y * sumy
		X2 = X2 + A2x * sumx + A2y * sumy
		j = j + 1
	end
	
	--
	--
	--

	local detc1c2 = C11 * C22 - C12 * C12
	local alpha_l, alpha_r = 0, 0
	
	if 1 + detc1c2 * detc1c2 ~= 1 then
		alpha_l = (C22 * X1 - C12 * X2) / detc1c2
		alpha_r = (C11 * X2 - C12 * X1) / detc1c2
	end

	--
	--
	--

	local seg_length = hypot(vl.x - vf.x, vl.y - vf.y)
	local epsilon = 1e-6 * seg_length
	
	if alpha_l < epsilon or alpha_r < epsilon then
		local dist = seg_length / 3

		return dist, dist
	else
		return alpha_l, alpha_r
	end
end

--
--
--

local Curve = { false, {}, {}, false }

--
--
--

local function AccumulateAndReturn (result, n)
	for i = 2, 4 do
		local cp = Curve[i]

		n, result[n + 1], result[n + 2] = n + 2, cp.x, cp.y
	end

	return n
end

--
--
--

local Param = {}

--
--
--

local function ChordLengthParam (points, first, last)
	Param[1] = 0
	
	local x1, y1, n, s = points[first].x, points[first].y, 1, 0
	
	for i = first + 1, last do
		local p = points[i]
		local x2, y2 = p.x, p.y

		s, n = s + hypot(x2 - x1, y2 - y1), n + 1

		Param[n], x1, y1 = s, x2, y2
	end

	for i = 2, n do
		Param[i] = Param[i] / s
	end
end

--
--
--

local function CenterTangent (points, index)
	local prev, next = points[index - 1], points[index + 1]

	return normalize(prev.x - next.x, prev.y - next.y)
end

--
--
--

local function OneSidedTangent (points, index, offset)
	local p1, p2 = points[index], points[index + offset]

	return normalize(p2.x - p1.x, p2.y - p1.y)
end

--
--
--

local function LeftTangent (points, index)
	return OneSidedTangent(points, index, 1)
end

--
--
--

local function RightTangent (points, index)
	return OneSidedTangent(points, index, -1)
end

--
--
--

local function ComputeMaxError (points, i1, i2, curve, u)
	local max_dist_sq, j, split = 0, 1
	
	for i = i1 + 1, i2 - 1 do
		local p, px, py = points[i], ComputeParam(3, curve, u[j])
		local vx, vy = px - p.x, py - p.y
		local dist_sq = vx * vx + vy * vy

		if dist_sq >= max_dist_sq then
			max_dist_sq, split = dist_sq, i
		end
		
		j = j + 1
	end
	
	return max_dist_sq, split or .5 * (i2 - i1 + 1) % 1
end

--
--
--

local function SetControlPoint (p, cp, tx, ty, alpha)
	cp.x, cp.y = p.x + tx * alpha, p.y + ty * alpha
end

--
--
--

local function FitCubic (points, i1, i2, tx1, ty1, tx2, ty2, err, result, n)
	local p1, p4 = points[i1], points[i2]
		
	Curve[1], Curve[4] = p1, p4

	if i2 == i1 + 1 then
		local dist = hypot(p4.x - p1.x, p4.y - p1.y) / 3

		SetControlPoint(p1, Curve[2], tx1, ty1, dist)
		SetControlPoint(p4, Curve[3], tx2, ty2, dist)

		return AccumulateAndReturn(result, n)
	end
	
	ChordLengthParam(points, i1, i2)

	local a1, a2 = Generate(points, i1, i2, Param, tx1, ty1, tx2, ty2)
	
	-- Find max deviation of points to fit curve.
	SetControlPoint(p1, Curve[2], tx1, ty1, a1)
	SetControlPoint(p4, Curve[3], tx2, ty2, a2)

	local max_error, split = ComputeMaxError(points, i1, i2, Curve, Param)
	
	if max_error < err then
		return AccumulateAndReturn(result, n)
	end
	
	-- If error not too large, try some re-parameterizations and iterations.
	if max_error < err ^ 2 then
		for _ = 1, 4 do
			for i = _, i2 - i1 + 1 do
				Param[i] = FindRoot(Curve, points, Param[i])
			end

			a1, a2 = Generate(points, i1, i2, Param, tx1, ty1, tx2, ty2)

			SetControlPoint(p1, Curve[2], tx1, ty1, a1)
			SetControlPoint(p4, Curve[3], tx2, ty2, a2)

			max_error, split = ComputeMaxError(points, i1, i2, Curve, Param)
			
			if max_error < err then
				return AccumulateAndReturn(result, n)
			end
		end
	end 	
	
	-- Fitting failed -> split, and recursion.
	local tcx, tcy = CenterTangent(points, split)

	n = FitCubic(points, i1, split, tx1, ty1, tcx, tcy, err, result, n)
	n = FitCubic(points, split, i2, -tcx, -tcy, tx2, ty2, err, result, n)

	return n
end

--
--
--

function M.FitCurve (points, n, err, result)
	local ltx, lty = LeftTangent(points, 1)
	local rtx, rty = RightTangent(points, n)

	result = result or {}
	result[1], result[2] = points[1].x, points[1].y

	return result, FitCubic(points, 1, n, ltx, lty, rtx, rty, err, result, 2)
end

--
--
--

return M