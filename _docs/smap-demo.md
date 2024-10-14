---
layout: default
title:  Tutorial
nav_order: 4

---

#  (Project Title) Tutorial

---

# (Example Below, please edit as per your project)

## AWS Deployment of R Shiny Bioinformatics App 
> An example of [sMAP: Standard Microarray Analysis Pipeline](https://github.com/BI-STEM-Away/sMAP).

### Launching an EC2 Instance from your AWS Management Console dashboard

A snapshot demo of launching an AWS EC2 Instance is described [here](AWS_EC2.pdf).

### Installing Docker in your EC2 Instance:

```Shell
sudo apt-get update

sudo apt-get install     ca-certificates     curl     gnupg     lsb-release

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

echo   "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com\
/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io
```


### Pulling and Running sMAP application in background at port 80 of your EC2 Instance:

`sudo docker run --rm -d -p 80:80 samuelbharti/smap:latest`


### Post-deployment: Demo of sMAP

A snapshot demo of sMAP analysis is described [here](sMAP_demo.pdf).

#### Workshop Hosted by STEM-AWAY
> Check out more of our internship and career advancement programs, workshops and mini-projects at [stemaway.com](https://stemaway.com/).


