-- premake5.lua
defines { "GLOBAL" }
	kind "StaticLib"
	language "C++"
	cppdialect "C++20"
	architecture "x86_64"
	location "build"

	objdir "%{cfg.location}/obj/%{cfg.platform}_%{cfg.buildcfg}"	
	targetdir "build/%{cfg.buildcfg}"

	includedirs { "../utils/src" }	
	includedirs { "../logger/src" }	
	includedirs { "../eilig/src" }	
	includedirs { "../basis/src" }
	includedirs { "../node/src" }	
	includedirs { "../kdtree/src" }		

workspace "interpolation"
	configurations { "Debug", "Release", "ReleaseCL" }

project "interpolation"

	files "src/interpolation.hpp"
	files "src/interpolation_types.hpp"
	files "src/interpolation_fixed.hpp"
	files "src/interpolation_fixed.cpp"
	files "src/interpolation_idw.hpp"
	files "src/interpolation_idw.cpp"
	files "src/interpolation_nn.hpp"
	files "src/interpolation_nn.cpp"
	files "src/interpolation_rbf.hpp"
	files "src/interpolation_rbf.cpp"

	filter "configurations:Debug"
		defines { "DEBUG" }
		symbols "On"

	filter "configurations:Release"   
		defines { "NDEBUG" }
		optimize "Speed"


	filter "configurations:ReleaseCL"   
		defines { "NDEBUG", "EILIG_ENABLE_OPENCL" }
		optimize "Speed"	
		
		
		files "src/interpolation_opencl_kernel.hpp"
		files "src/interpolation_opencl_kernel.cpp"
		files "src/interpolation_opencl_rbf.hpp"
		files "src/interpolation_opencl_rbf.cpp"
		
		includedirs { "../club/src" }	
		includedirs { "../opencl/inc" }