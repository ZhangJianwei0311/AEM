all : AEM

AEM : AEMtrain.o topiclib.o
	g++ -g -o AEM AEMtrain.o topiclib.o

AEMtrain.o : AEMtrain.cpp topiclib.h
	g++ -g -c AEMtrain.cpp

topiclib.o : topiclib.cpp topiclib.h
	g++ -g -c topiclib.cpp

clean :
	rm AEM AEMtrain.o topiclib.o
