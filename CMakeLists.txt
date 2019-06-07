cmake_minimum_required(VERSION 3.8)
project(density CXX)
# 最適化・警告等のオプション
set(CMAKE_CXX_STANDARD 17)

# curdirをinclude file探索対象にする
# mocで生成されたファイルをincludeするはずなので、有効にしとくが吉
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# mocコマンドを自動実行する
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Qt5ライブラリを検索する。
find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)

# Qt5のインストール先フォルダをcmake_prefix_pathに追加(find_packageの探索先)
# cmakeのmodule配下にqtがインストールされてる場合は必要なし
set(Qt5_DIR /usr/local/Cellar/qt)
list(APPEND CMAKE_PREFIX_PATH ${Qt5_DIR})

file(GLOB_RECURSE CPP_FILES *.cpp)

include_directories(.)

add_executable(density density.cpp common.hpp plot.hpp plot.cpp load_json.cpp)
add_executable(biterror biterror.cpp common.hpp plot.hpp plot.cpp load_json.cpp)

set(CMAKE_CXX_FLAGS "-O3 -Wall -g -fPIC -fmessage-length=0 -I/usr/local/include -I/usr/local/Cellar/qt/5.12.3/include -I/usr/local/Cellar/qwt/6.1.3_4/lib/qwt.framework/Versions/6/Headers -L/usr/local/lib")
target_link_libraries(density LINK_PUBLIC Qt5::Widgets Qt5::Core qwt fftw3 pthread m)
target_link_libraries(biterror LINK_PUBLIC Qt5::Widgets Qt5::Core qwt fftw3 pthread m)