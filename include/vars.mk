OBJDIR = mod
BINDIR = bin
SRCDIR = src
TESTDIR = test
SRCDIRS = src/assemble src/convert src/general src/misc src/orth src/orth_rank \
	src/qr_factorization src/row_compress src/solve src/transforms src/types
FC = gfortran
LFLAGS =
CFLAGS = -fbounds-check -mcmodel=medium -Waliasing -Wampersand -Wconversion -Wsurprising \
	-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wunused
#CFLAGS = -Wall -O2 -mcmodel=medium
PROFLFLAGS = -pg
PROFCFLAGS = -pg
