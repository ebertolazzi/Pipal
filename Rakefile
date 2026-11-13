# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Maximum number of parallel compile jobs. GCC/MSVC can need a couple of GB per
# C++20 translation unit, and an over-subscribed build is a classic source of
# compiler crashes. PINS_MAX_CORE (or JOBS) caps it; 16 by default.
def pins_jobs
  cap = ENV['PINS_MAX_CORE'].to_i
  cap = 16 if cap <= 0
  n = ENV['JOBS'].to_i
  n = Etc.nprocessors if n <= 0
  [n, cap].min
end

# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Pipal project is distributed under the MIT License.                                         #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

# Configuration of the build
BUILD_DEBUG = false
BUILD_TESTS = true

case RUBY_PLATFORM
when /mingw|mswin/
  PARALLEL = ""
  QUIET    = ""
else
  require 'etc'
  cmakeversion = %x( cmake --version ).scan(/\d+\.\d+\.\d+/).last
  mm = cmakeversion.split('.');
  if mm[0].to_i > 3 || (mm[0].to_i == 3 && mm[1].to_i >= 12) then
    PARALLEL = "--parallel #{pins_jobs} "
    QUIET    = ""
  else
    PARALLEL = ""
    QUIET    = ""
  end
end

# MSYSTEM is exported by every MSYS2 shell and is the only reliable way to
# tell a MinGW environment from a plain POSIX one: the ruby shipped in the
# MSYS subsystem reports its platform as "linux"/"cygwin" even when it drives
# a UCRT64 / MINGW64 toolchain, which made this build announce itself as Linux
# and pick the Linux library suffixes. MSYSTEM=MSYS is the POSIX shell and is
# deliberately NOT treated as MinGW.
PINS_MINGW_ENVS = %w(UCRT64 MINGW64 MINGW32 CLANG64 CLANG32 CLANGARM64)
PINS_IS_MINGW   = PINS_MINGW_ENVS.include?(ENV['MSYSTEM'].to_s.upcase)

# The previous dispatch was inverted: a Windows-looking RUBY_PLATFORM selected
# :build_linux and everything else fell through to :build_windows, so a MSYS2
# ruby (which looks like linux) built with the Windows recipe.
if PINS_IS_MINGW || (/mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
  task :default => [:build_windows]
elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
  task :default => [:build_osx]
else
  task :default => [:build_linux]
end

cmd_cmake_build = "-G Ninja "
if BUILD_TESTS then
  cmd_cmake_build += "-DPIPAL_BUILD_TESTS:VAR=true "
else
  cmd_cmake_build += "-DPIPAL_BUILD_TESTS:VAR=false "
end
if BUILD_DEBUG then
  cmd_cmake_build += "-DCMAKE_BUILD_TYPE:VAR=Debug "
else
  cmd_cmake_build += "-DCMAKE_BUILD_TYPE:VAR=Release "
end

task :default => [:build]

TESTS = []

desc "Run tests"
task :run do
  puts "run test".yellow
  Dir.glob('./build/tests/test_*') do |cmd|
    next if cmd =~ /.manifest$|.dSYM$/
    puts "execute: #{cmd}".yellow
    sh cmd
  end
end

task :build_gen do

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build

  puts "run CMake for Pipal".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMake for Pipal".yellow
  if BUILD_DEBUG then
    sh 'cmake --build . --config Debug ' + PARALLEL + QUIET
  else
    sh 'cmake --build . --config Release ' + PARALLEL + QUIET
  end
  FileUtils.cd '..'
end

desc "Compile for OsX"
task :build_osx => :build_gen do |t, args|
end

desc "Compile for Linux"
task :build_linux => :build_gen do |t, args|
end

desc "Compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build

  puts "run CMake for Pipal".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMake for Pipal".yellow
  if BUILD_DEBUG then
    sh 'cmake --build . --config Debug ' + PARALLEL + QUIET
  else
    sh 'cmake  --build . --config Release  ' + PARALLEL + QUIET
  end

  FileUtils.cd '..'
end

desc "Build for OsX/Linux/Windows"
task :build do
  if PINS_IS_MINGW || (/mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
    Rake::Task["build_windows"].invoke
  elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
    Rake::Task["build_osx"].invoke
  else
    Rake::Task["build_linux"].invoke
  end
end

task :clean_gen do
  FileUtils.rm_rf 'bin'
  FileUtils.rm_rf 'build'
  FileUtils.rm_rf 'third_party'
end

desc "Clean for OsX"
task :clean_osx => :clean_gen do
end

desc "Clean for Linux"
task :clean_linux => :clean_gen do
end

desc "Clean for Windows"
task :clean_win => :clean_gen do
  FileUtils.rm_rf 'vs_*'
end

desc "Clean for OsX/Linux/Windows"
task :clean do
  if (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
    # Linux
    Rake::Task["clean_linux"].invoke
  elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
    # OsX
    Rake::Task["clean_osx"].invoke
  else
    # Windows
    Rake::Task["clean_win"].invoke
  end
end
