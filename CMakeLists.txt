cmake_minimum_required(VERSION 3.8)
project(density CXX)
add_executable(density density.cpp check_node.cpp variable_node.cpp fourier_gauss.hpp transform_complex.hpp common.hpp QWindows.hpp)
# 最適化・警告等のオプション
set(CMAKE_CXX_STANDARD 17)

# curdirをinclude file探索対象にする
# mocで生成されたファイルをincludeするはずなので、有効にしとくが吉
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# mocコマンドを自動実行する
set(CMAKE_AUTOMOC ON)

# Qt5のインストール先フォルダをcmake_prefix_pathに追加(find_packageの探索先)
# cmakeのmodule配下にqtがインストールされてる場合は必要なし
set(Qt5_DIR /usr/local/Cellar/qt)
list(APPEND CMAKE_PREFIX_PATH ${Qt5_DIR})

# 適当に使うモジュールを選択
# Qt5<module name>
find_package(Qt5Widgets REQUIRED)

set(CMAKE_CXX_FLAGS "-O3 -Wall -g -fPIC -fmessage-length=0 -I/usr/local/include -I/usr/local/Cellar/qt/5.12.2/include -L/usr/local/lib")
target_link_libraries(density LINK_PUBLIC Qt5::Widgets fftw3 pthread m)