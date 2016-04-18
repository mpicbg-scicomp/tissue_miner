System requirement
================

* **Linux**

Docker requires a 64-bit installation regardless of your Linux version and distribution. Additionally, your kernel must be 3.10 at minimum. The latest 3.10 minor version or a newer maintained version are also acceptable.
See installation procedure for ubuntu [here](https://docs.docker.com/engine/installation/ubuntulinux/)

* **MacOSX**

Your Mac must be running OS X 10.8 “Mountain Lion” or newer to install the Docker Toolbox.
See installation procedure [here](https://docs.docker.com/engine/installation/mac/)

* **Windows**

Your machine must be running Windows 7 or newer to run Docker. 
See installation procedure [here](https://docs.docker.com/engine/installation/windows/)


"Cannot connect to the Docker daemon. Is the docker daemon running on this host?"
================

This may happen on older MacOS installation prior to Yosemite.
We recommend to follow the following procedure:

     ##  In the shell launched by Docker Quickstart Terminal, list available virtual machines 
     docker-machines ls
    
     ## If the list is empty then create a default else go to the next step  
     docker-machine create --driver virtualbox default
    
     ## Connect to the virtual machine called default
     docker-machine ssh default
    
     ## Pull the TissueMiner application image
     docker pull brandl/tissue_miner
    
     ## Go to your movie repository  and run the TissueMiner automated workflow on the Demo 

[see How to get started](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/README.md#how-to-get-started)
    
