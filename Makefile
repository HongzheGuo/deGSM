
all : deGSM ubwt
.PHONY : all

CC = gcc
CFLAGS = -g -Wall -O2 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
#CFLAGS = -g -Wall -O0 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
LIB = -lz -lpthread -lm

#BIN_DIR = .
SRC_DIR_deGSM = ./src_deGSM

SOURCE_deGSM = $(wildcard ${SRC_DIR_deGSM}/*.c) 

OBJS_deGSM = $(SOURCE_deGSM:.c=.o)

SRC_DIR_ubwt = ./src_ubwt

SOURCE_ubwt = $(wildcard ${SRC_DIR_ubwt}/*.c) 

OBJS_ubwt = $(SOURCE_ubwt:.c=.o)

deGSM : $(OBJS_deGSM) 
	$(CC) $(CFLAGS) -o deGSM $(OBJS_deGSM) $(LIB)

ubwt : $(OBJS_ubwt) 
	$(CC) $(CFLAGS) -o ubwt $(OBJS_ubwt) $(LIB)	
	
.PHONY:clean
clean : 
	rm deGSM $(OBJS_deGSM) ubwt $(OBJS_ubwt) 
