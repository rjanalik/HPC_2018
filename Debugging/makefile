CFLAGS = 
#CFLAGS = -g
#CFLAGS = -g -fsanitize=address

all: dotProduct

dotProduct: dotProduct.cpp
	g++ $(CFLAGS) $< -o $@

clean:
	rm -rf dotProduct
