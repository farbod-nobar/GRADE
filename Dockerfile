FROM alpine:latest

# install GCC 6.10 and other required packages
RUN apk update && apk add gcc g++ make vim

# Copy and build your C++ code
COPY . /Grade
WORKDIR /Grade
RUN ["make"]
RUN ["cp", "GRADE", "/usr/local/bin/grade"]



