--- TODO

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

-- Exports --
local M = {}

--
--
--

do
	local Row = {}

	--- Subdivide a B&eacute;zier spline into two new splines, using [De Casteljau's algorithm](http://en.wikipedia.org/wiki/De_Casteljau's_algorithm).
	-- @array bezier **Vector** elements 1, 2, ..., _deg_ + 1 corresponding to the first
	-- endpoint, first control point, ..., final endpoint of the B&eacute;zier spline to subdivide.
	--
	-- It is safe to reuse the value of _bezier_ as either _dst1_ or _dst2_.
	-- @array dst1 Receives the "left" subdivision, i.e. the spline is evaluated from 0 to
	-- _t_ and **Vector** elements 1, 2, ..., _deg_ are populated with the results (tables
	-- being created if necessary).
	-- @array dst2 As per _dst1_, the "right" subdivision, evaluated from _t_ to 1.
	-- @number[opt=.5] t Parameter at which to split the spline.
	-- @int[opt=#bezier - 1] deg Degree of the spline.
	function M.Subdivide (bezier, dst1, dst2, t, deg)
		t = t or .5

		-- The base of the De Casteljau triangle is just the source spline.
		local n, height = 0, deg and deg + 1 or #bezier

		for i = 1, height do
			Row[n + 1] = bezier[i].x
			Row[n + 2] = bezier[i].y

			n = n + 2
		end

		-- Iterate up the triangle. The left-hand side, top to bottom, supplies the left
		-- subdivision; the right-hand side, bottom to top, supplies the right.
		local s = 1 - t

		for i = 1, height do
			local j = height - i + 1
			local vl = dst1[i] or {}
			local vr = dst2[j] or {}

			vl.x, vr.x = Row[1], Row[n - 1]
			vl.y, vr.y = Row[2], Row[n - 0]

			dst1[i], dst2[j] = vl, vr

			-- Generate the next row up by interpolating the working row. This can be
			-- performed in place, since the left-hand interpoland will not be needed
			-- after the current step and thus can be overwritten.
			for k = 1, n - 2, 2 do
				Row[k + 0] = s * Row[k + 0] + t * Row[k + 2]
				Row[k + 1] = s * Row[k + 1] + t * Row[k + 3]
			end

			n = n - 2
		end
	end
end

--
--
--

return M