#
# Makefile
#
CC      = g++
CFLAGS  = -O3
PROG    = AdapterRemoval


all: $(PROG) man

$(PROG): $(PROG).cc
	$(CC) $(PROG).cc $(CFLAGS) -o $(PROG)

man: 
	pod2man $(PROG).pod > $(PROG).1

# Clean
clean:
	rm -f $(PROG) $(PROG).1 *~

# Install
install:
	mv -f $(PROG) /usr/local/bin/
	mv -f $(PROG).1 /usr/share/man/man1/
