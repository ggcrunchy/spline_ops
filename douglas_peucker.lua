--- Adapted from **quangnle**'s implementation of the Douglas-Peucker algorithm, as found in [CurveFitting](https://github.com/quangnle/CurveFitting_PJSchneider).
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

-- Standard library imports --
local abs = math.abs

-- Exports --
local M = {}

--
--
--

local function AuxReduce (points, i1, i2, tolerance, indices, ni)
	if i2 - i1 > 1 then -- not just endpoints?
		local p1, p2 = points[i1], points[i2]
		local px, py = p1.x, p1.x
		local wx, wy = p2.x - px, p2.y - py
		local denom = wx * wx + wy * wy

		if 1 + denom ~= 1 then -- segment not empty? (TODO: just the loop case at a lower level?)
			local wxp, wyp, max_dist_sq, imax = wx / denom, wy / denom, 0
			
			--
			--
			--
			
			for i = i1 + 1, i2 - 1 do -- ignore endpoints
				local p = points[i]
				local vx, vy = p.x - px, p.y - py
				local t = vx * wxp + vy * wyp
				local rx, ry = vx - t * wx, vy - t * wy
				local dist_sq = rx * rx + ry * ry
			
				if dist_sq > max_dist_sq then
					max_dist_sq, imax = dist_sq, i
				end
			end
			
			--
			--
			--
			
			if max_dist_sq > tolerance and imax ~= 1 and imax ~= i2 then
				indices[ni + 1], ni = imax, ni + 1

				ni = AuxReduce(points, i1, imax, tolerance, indices, ni)
				ni = AuxReduce(points, imax, i2, tolerance, indices, ni)
			end
		end
	end

	return ni
end

--
--
--

local Indices = {}

--
--
--

--- DOCME
function M.Reduce (points, n, tolerance, out)
	if not points or n < 3 then
		return points, points and n or 0
	end

	--
	--
	--
	
	local p1, plast = points[1], points[n]
	local pn, p1x, p1y = plast, p1.x, p1.y 
	
	while pn ~= p1 and abs(p1x - pn.x) < 1e-3 and abs(p1y - pn.y) < 1e-3 do
		n = n - 1
		pn = points[n]
	end

	local was_loop = pn ~= plast

	--
	--
	--

	Indices[1], Indices[2] = 1, n

	--
	--
	--

	local ni = AuxReduce(points, 1, n, tolerance, Indices, 2)
	
	--
	--
	--
	
	for i = 1, ni - 1 do
		for j = ni, i + 1, -1 do
			if Indices[i] > Indices[j] then
				Indices[i], Indices[j] = Indices[j], Indices[i]
			end
		end
	end
	
	--
	--
	--
	
	out = out or {}
	
	for i = 1, ni do
		out[i] = points[Indices[i]]
	end
	
	return out, ni, was_loop
end

--
--
--

return M