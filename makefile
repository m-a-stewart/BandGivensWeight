include ./include/vars.mk
OBJDIR = mod
test : force_look
	cd test; $(MAKE)

force_look :
	true

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TEST_EXECS)
