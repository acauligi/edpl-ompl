#!/bin/bash

# update gcc/g++
sudo apt install -y gcc-5 g++-5
sudo apt install -y git cmake libeigen3-dev mercurial 
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5

# install IDE
# sudo apt install codeblocks

# set up the source code directory
mkdir -p ~/ros/firm_ws/src
cd ~/ros/firm_ws/src

# install source code
git clone https://github.com/sauravag/edpl-ompl.git

# install dependencies
	# install omplapp from source
	git clone https://github.com/ompl/omplapp.git
  cd omplapp
	git checkout origin/1.3.1	# d0511
		# install ompl from source
			git clone https://github.com/ompl/ompl.git
      cd ompl # acauligi 
			git checkout origin/1.3.1	# ab02d
      cd .. # acauligi
				# install_common_dependencies()
					sudo apt install -y libboost-all-dev libflann-dev
					sudo apt install -y libode-dev
					sudo apt install -y python-pip
				# install_python_binding_dependencies()
					sudo -H pip install -vU pygccxml https://bitbucket.org/ompl/pyplusplus/get/1.8.0.tar.gz
					#sudo apt install castxml
						cd ~
						# wget https://midas3.kitware.com/midas/download/item/318227/castxml-linux.tar.gz
            wget https://data.kitware.com/api/v1/file/5877a2498d777f05f44b13fb/download/castxml-linux.tar.gz   # acauligi
						sudo apt install stow	# in case you don't have this
						sudo tar zxf castxml-linux.tar.gz -C /usr/local/stow
						rm castxml-linux.tar.gz
						cd -
						cd /usr/local/stow
						sudo stow castxml
						cd -
				# build ompl
					mkdir -p ompl/build/Release
					cd ompl/build/Release
					cmake ../..
					#make update_bindings
					make -j8
					sudo make install
					cd -
		# install_app_dependencies()
		cd ~/ros/firm_ws/src/omplapp
			ubuntu_version=`lsb_release -rs | sed 's/\.//'`
			mkdir omplapp_dep
			cd omplapp_dep
			# We prefer PyQt5, but PyQt4 also still works.
				if [[ $ubuntu_version > 1410 ]]; then
					sudo apt install -y python-pyqt5.qtopengl
				else                 
					sudo apt install -y python-qt4-dev python-qt4-gl
				fi
				sudo apt install -y freeglut3-dev libassimp-dev python-opengl python-flask python-celery
			# install additional python dependencies via pip
				sudo -H pip install -vU PyOpenGL-accelerate
			# install libccd
				if [[ $ubuntu_version > 1410 ]]; then
					sudo apt install -y libccd-dev
				else
					wget -O - https://github.com/danfis/libccd/archive/v2.0.tar.gz | tar zxf -
					rm v2.0.tar.gz
					cd libccd-2.0; cmake .; sudo -E make install; cd ..
				fi
			# install fcl from source
				git clone https://github.com/flexible-collision-library/fcl.git
				cd fcl
				git checkout 0.3.2
				mkdir -p build
				cd build
				cmake ..
				make -j8
				sudo make install
				cd ../..
			# install fcl from source
# 				if ! pkg-config --atleast-version=0.5.0 fcl; then
# 					if [[ $ubuntu_version > 1604 ]]; then
# 						sudo apt install -y libfcl-dev
# 					else
# 						wget -O - https://github.com/flexible-collision-library/fcl/archive/0.5.0.tar.gz | tar zxf -
# 						rm 0.5.0.tar.gz
# 						cd fcl-0.5.0; cmake .; sudo -E make install; cd ..
# 					fi
# 				fi
			cd ..
		# build omplapp
			mkdir -p build/Release
			cd build/Release
			cmake ../..
			#make update_bindings
			make -j8
			sudo make install
			cd -
		cd ..

	# install armadillo from debian package
# 	sudo apt-get install libarmadillo-dev
	# install armadillo from source
		sudo apt install -y libopenblas-dev liblapack-dev libarpack2-dev
		# wget http://sourceforge.net/projects/arma/files/armadillo-7.900.1.tar.xz
    wget https://sourceforge.net/projects/arma/files/armadillo-7.950.0.tar.xz # acauligi
		tar -xvf armadillo-7.*
		rm armadillo-7.*.tar.xz
		mv armadillo-7.* armadillo
		cd armadillo
		cmake . -DCMAKE_INSTALL_PREFIX=/usr/local/stow/armadillo
		make -j8 
		sudo make install
		cd -
		cd /usr/local/stow
		sudo stow armadillo
		cd -

# install other dependencies
	sudo apt install -y freeglut3-dev libqt4-dev libtinyxml-dev
	#sudo apt-get install ros-indigo-aruco-ros

# build firm package
cd edpl-ompl
git checkout firmcp-dev
#git checkout baseline
# git checkout 8e54b	# this is the stable version for URM-POMCP
mkdir -p build && cd build
cmake ..
make -j8
