# Graph-Mp
## Introduction
Graph-MP (A Generalized Matching Pursuit Approach for Graph-Structured Sparsity). We already cited the authors of all related open source code. If you use our code, please cite our paper "A Generalized Matching Pursuit Approach for Graph-Structured Sparsity" (<a href="http://www.ijcai.org/Proceedings/16/Papers/200.pdf">link</a>), which will be appeared in IJCAI-2016. I want to say thank you to Ludwig Schmidt who is a Ph.D. student from MIT. He provided his fast-pcst c++ code to me. You can access the fast-pcst c++ code in the this <a href="https://github.com/ludwigschmidt/pcst-fast">link</a>

## Before you go
Before you run this code, make sure you have maven and Java 1.8 on your machine. I tested all of the code under Linux environment. If you have errors or problems when you use the code on Windows or Mac OS, please let me know. My email is bzhou6@albany.edu

## Dependency Information
Need to be installed:
#### 1. <a href="https://maven.apache.org/install.html" >Apache Maven 3.3.x</a>
#### 2. <a href="http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html">Java 1.8</a>
#### 3. <a href="https://git-scm.com/book/en/v2/Getting-Started-Installing-Git">git</a>

## How to run the code ?
#### Step-1: download the code using the following command:
```shell
git clone https://github.com/baojianzhou/Graph-MP.git
```
#### Step-2: go into the folder path/to/Graph-MP.
```shell
cd path/to/Graph-MP
```
#### Step-3: compile the code under Java 1.8
```shell
mvn compile
```
#### Step-4: test sample testing cases
```shell
mvn test
```
