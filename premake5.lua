-- premake5.lua
workspace "interpolation"
	configurations { "Debug", "Release" }
	location "build"

project "interpolation"
	kind "StaticLib"
	language "C++"
	cppdialect "C++17"

	targetdir "build/%{cfg.buildcfg}"
	includedirs { "../utils/src" }	
	includedirs { "../logger/src" }	
	includedirs { "../eilig/src" }	
	includedirs { "../basis/src" }
	includedirs { "../node/src" }	
	includedirs { "../kdtree/src" }	

	files { "src/**.hpp", "src/**.cpp" }

	filter "configurations:Debug"
		architecture "x86_64"     
		defines { "DEBUG" }
		symbols "On"

	filter "configurations:Release"
		architecture "x86_64"     
		defines { "NDEBUG" }
		optimize "Speed"