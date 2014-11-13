# encoding: utf-8
# -*- mode: ruby -*-
# vi: set ft=ruby :
# This VagrantFile is used to setup a guest virtual machine on your computer.
# type:
#    $> vagrant up
# to start the vm. (get vagrant from http://www.vagrantup.com)

# Vagrant.require_version ">= 1.2.2"
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "MotionDetection-VM"
  config.vm.box_url = "http://cloud-images.ubuntu.com/vagrant/trusty/current/"\
                      "trusty-server-cloudimg-amd64-vagrant-disk1.box"
  config.ssh.forward_agent = true

  # Map through port 8888 for IPython Notebook
  config.vm.network "forwarded_port", guest: 8888, host: 8888

  # install required packages in guest
  config.vm.provision "shell", inline: "apt-get update"
  config.vm.provision "shell", inline: "apt-get install --yes"\
            " texlive-fonts-recommended texlive-latex-extra"\
            " texlive-latex-base python-progressbar ipython"\
            " python-matplotlib ipython-notebook python-numpy"\
            " python-scipy python-sympy supervisor inkscape"\
            " pandoc jq python-h5py libzmq3-dev python-pip"
  config.vm.provision "shell", inline: "pip install runipy"

  config.vm.provision "shell", inline: "cd /vagrant/tutorial-src && make html html-basic"

  # automatically startup ipython 
  config.vm.provision "shell", inline: "cat >"\
            "/etc/supervisor/conf.d/ipython.conf <<EOF\n"\
            "[program:ipython-notebook]\n"\
            "# This will run an ipython web notebook at http://0.0.0.0:8888\n"\
            "command=ipython notebook --NotebookApp.ip='*' --no-browser --NotebookApp.port=8888 --matplotlib=inline\n"\
            "environment=HOME='/vagrant/tutorial-src', USER='vagrant'\n"\
            "directory=/vagrant/tutorial-src/\n"\
            "EOF\n"\
            "supervisorctl reload"

  # Once this is done running, you can open the IPython notebook at
  # http://127.0.0.1:8888/ .

end
