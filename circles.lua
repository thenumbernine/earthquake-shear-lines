local vec3d = require 'vec-ffi.vec3d'

local M = {}

local dayInSec = 60 * 60 * 24


-- when searching for great-arc circles, ignore point-pairs that are not within this angle distance apart
--M.greatArcAngleMin, M.greatArcAngleMax = 10, 90
--M.greatArcAngleMin, M.greatArcAngleMax = 0, 180
M.greatArcAngleMin, M.greatArcAngleMax = 10, 170

-- tolerance of how close a quake must be to a previous circle to consider along the geodesic
M.quakeAlignWithPreviousGeodesicAngleThreshold = 1	-- in radians
--M.quakeAlignWithPreviousGeodesicAngleThreshold = 5	-- in radians
--M.quakeAlignWithPreviousGeodesicAngleThreshold = 15	-- in radians

-- if we get two great-arcs touching our current quake (within `quakeAlignWithPreviousGeodesicAngleThreshold`)
--  then throw out the new one if their axii (+-) is within this angle apart.
-- ... but I think we do want duplicate-arcs to show where multiple-quakes-on-same-arc do align.
-- ... we just don't want duplicate-points, being used to generate those arcs... right?
M.filterDuplicateTouchingArcsAngleThreshold = 0
--M.filterDuplicateTouchingArcsAngleThreshold = 1
--M.filterDuplicateTouchingArcsAngleThreshold = 5

-- when searching for great-arc circles, ignore points within this many radians of previously-considered points
--M.filterDuplicatePointsAngleThreshold = 1
--M.filterDuplicatePointsAngleThreshold = 5
M.filterDuplicatePointsAngleThreshold = 10

-- how far back we want to look when comparing great-arcs 
--M.timeWindowInDays = 1
M.timeWindowInDays = 3


local basisFor = |v| do
	local x = v:cross(vec3d(1, 0, 0))
	local y = v:cross(vec3d(0, 1, 0))
	local z = v:cross(vec3d(0, 0, 1))
	local xl = x:lenSq()
	local yl = y:lenSq()
	local zl = z:lenSq()
	-- TODO there's too many situations where ternary needs ()'s to save parsing... I should just change its symbols
	return xl > yl	-- x > y
		and (
			xl > zl
			and x	-- x > y, x > z
			or z		-- z > x > y
		)
		or (		-- y >= x
			yl > zl
			and y	-- y > z, y >= x
			or z	-- z > y >= x
		)
end

M.calcCircles=|earthquakes|do
	local activeCircles = table()	-- holds indexes into 'allCircles'
	local allCircles = table()
	local oldestEarthquakeIndex = 1
	for i,eq in ipairs(earthquakes) do
		local v = eq.xyz
		
		-- see how close it is to all other circles on file
		local touchingCircles = table()
		for ici,ci in ipairs(activeCircles) do
			local c = allCircles[ci]
			local angle = math.acos(c.axis:dot(v))
			-- if we were right on the plane then we would be right at 90' with the axis.
			-- see how far the angle we make with the axis is from 90'
			if math.abs(angle - .5 * math.pi) < math.rad(M.quakeAlignWithPreviousGeodesicAngleThreshold) then
				
				-- if a previous great-arc was already associated with this point and it is close enough to the current great-arc then skip the current one (uniques only)
				local dups
				if M.filterDuplicateTouchingArcsAngleThreshold then
					local costh = math.cos(math.rad(M.filterDuplicateTouchingArcsAngleThreshold))
					for ici2=1,ici-1 do
						local c2 = allCircles[activeCircles[ici2]]
						if math.abs(c.axis:dot(c2.axis)) > costh then
							dups = true
							break
						end
					end
				end
				if not dups then
					-- TODO take note (or even filter?) based on the angle distance between us and the prevoius points that it took to form this geodesic?
					-- what if all new quakes happen to be along geodesics but >90' of old geodesic points? (or whatever our geodesic creation angle criteria is)
					touchingCircles:insert{circle=ci, angle=angle}
				end
			end
		end
		local numCirclesChecked = #activeCircles

		-- now do our typical building of great-arcs based on our other angle criteria
		local oldCircleCount = #activeCircles
		for j=oldestEarthquakeIndex,i do
			-- TODO filter out repetitive points

			local eq2 = earthquakes[j]
			local v2 = eq2.xyz	-- spherical
			
			local dups 
			-- filter out duplicate points for great-arc construction consideration ...
			if M.filterDuplicatePointsAngleThreshold then
				local cosAngle = math.cos(math.rad(M.filterDuplicatePointsAngleThreshold))
				for k=oldestEarthquakeIndex,j-1 do
					local v3 = earthquakes[k].xyz	 -- spherical
					-- secant length to angle ...
					if v3:dot(v2) >= cosAngle then
						dups = true
						break
					end
				end
			end
			if not dups then
				local torque = v:cross(v2)	-- TODO times pair of earthquake magnitudes or something?
				local axisLen = torque:length()
				local axis = torque * (1 / (math.max(axisLen, 1e-15)))	
				local influence = axisLen * 10^eq.mag * 10^eq2.mag
				local angle = math.acos(v:dot(v2))
				if math.rad(M.greatArcAngleMin) < angle and angle < math.rad(M.greatArcAngleMax) then	-- if I weight by cross then that'll make the 0' and 180' angles diminish ... nah, there's still a lot of noise unless I turn down the alpha ... then it's hard to highlight any geodesics ...
					local a1 = axis
					local a2 = basisFor(a1):normalize()
					local a3 = a1:cross(a2)			
					local circle = {
						axis = axis,
						axis2 = a2,
						axis3 = a3,
						angle = math.pi * .5,
						indexes = {i, j},
						-- [[ used for rendering
						influence = influence,
						color = vec3d(1, .3, .07),
						centerPos = (v + v2):normalize(),	-- midpoint between two earthquake lines
						centerAngle = .5 * angle,			-- half the angle between them
						--]]
					}

					allCircles:insert(circle)
					activeCircles:insert(#allCircles)
				end
			end
		end
		-- now recalculate the oldest point index t consider (basd on our time window)
		while oldestEarthquakeIndex < i
		and eq.ostime - earthquakes[oldestEarthquakeIndex].ostime > M.timeWindowInDays * dayInSec
		do
			oldestEarthquakeIndex+=1 
		end
		-- and throw out old circles
		for j=oldCircleCount,1,-1 do
			local ci = activeCircles[j]
			local c = allCircles[ci]
			if c.indexes[1] < oldestEarthquakeIndex
			or c.indexes[2] < oldestEarthquakeIndex
			then
				activeCircles:remove(j)
			end
		end

		-- collect unique points to all circles of this point
		local uniquePointIndexes = {}
		uniquePointIndexes[i] = true
		for _,t in ipairs(touchingCircles) do
			local ci = t.circle
			local c = allCircles[ci]
			uniquePointIndexes[c.indexes[1]] = true
			uniquePointIndexes[c.indexes[2]] = true
		end
		eq.touchingCircles = touchingCircles
		eq.uniquePoints = uniquePoints	-- do I need a list of all points that are on all circles that are touching this point? or nah?

		--[[
		print('i='..i
			..' t='..eq.time
			..' touch='..#touchingCircles
			..' lat='..eq.latitude
			..' lon='..eq.longitude
			..' xyz='..eq.xyz
			..' mag='..eq.mag
			..' place='..eq.place
			..' check='..numCirclesChecked
			..' numPtsNow='..(i-oldestEarthquakeIndex+1)
			..' numArcsNow='..#activeCircles
		)
		--]]
		--[[ show what it's touching?
		print'\tcircles:'
		for _,t in ipairs(touchingCircles) do
			local ci = t.circle
			local c = allCircles[ci]
			print('\t\t'
				..' index='..c.indexes[1]..','..c.indexes[2]
				--..' index='..table.concat(c.indexes,',')	-- why would this segfault?
				..' angle='..t.angle
				..' axis='..c.axis
			)
		end
		print'\tpoints:'
		for _,j in ipairs(table.keys(uniquePointIndexes):sort()) do
			print('\t\t'..j..' '..earthquakes[j].xyz)
		end
		--]]
	end

	--[[ big and slow
	-- map out the vec3d's
	for _,eq in ipairs(earthquakes) do
		eq.xyz = {eq.xyz:unpack()}
	end
	for _,c in ipairs(allCircles) do
		c.axis = {c.axis:unpack()}
	end
	local luapath = path(basefn..'.lua')
	luapath:write(tolua{
		earthquakes = earthquakes,
		circles = allCircles,
	})
	--]]

	return allCircles
end

return M
