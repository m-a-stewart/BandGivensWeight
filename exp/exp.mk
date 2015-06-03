_EXPS = exp
EXPS = $(patsubst %,$(BINDIR)/%,$(_EXPS))

.PHONY : all_exps
all_exps : $(EXPS)

.PHONY : run_all_exps
run_all_exps : $(EXPS)
	$(BINDIR)/exp

.PHONY : run_exp
run_exp : $(BINDIR)/exp
	$(BINDIR)/exp

$(OBJDIR)/exp.o : $(OBJDIR)/orb.o $(OBJDIR)/test_data.o
