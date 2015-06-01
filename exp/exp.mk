_EXPS = exp1
EXPS = $(patsubst %,$(BINDIR)/%,$(_EXPS))

.PHONY : all_exps
all_exps : $(EXPS)

.PHONY : run_all_exps
run_all_exps : $(EXPS)
	$(BINDIR)/exp1

.PHONY : run_general
run_exp1 : $(BINDIR)/exp1
	$(BINDIR)/exp1

$(OBJDIR)/exp1.o : $(OBJDIR)/orb.o $(OBJDIR)/test_data.o
