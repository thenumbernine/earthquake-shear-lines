-- [=======[ TODO cache calcs somehow
local allCharts = require 'geographic-charts'
local chartNames = table{	-- which charts we want to allow ...
	'Equirectangular',
	'Azimuthal equidistant',
	'Mollweide',
	'WGS84',
	'sphere',
}
local charts = table()
for i,name in ipairs(chartNames) do
	local chart = allCharts![name]
	charts[i] = chart
	charts[name] = chart
end
for _,name in ipairs(chartNames) do
	charts[name]?:build()
end

local chartCNames = charts:mapi(|chart| chart:getCName())
local chartCode = require 'geographic-charts.code'(charts)
--]=======]

return {
	charts = charts,
	chartCode = chartCode, 
	chartCNames = chartCNames,
}
