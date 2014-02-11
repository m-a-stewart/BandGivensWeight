include ./include/vars.mk
OBJDIR = mod
test : force_look
	cd test; $(MAKE)

force_look :
	true

clean :
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/test_decomp $(OBJDIR)/test_convert \
		$(OBJDIR)/test_compress $(OBJDIR)/test_qr_iteration
