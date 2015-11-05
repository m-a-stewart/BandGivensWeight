_EXPS = exp exp1
EXPS = $(patsubst %,$(BINDIR)/%,$(_EXPS))

.PHONY : all_exps
all_exps : $(EXPS)

.PHONY : run_all_exps
run_all_exps : $(EXPS)
	$(BINDIR)/exp
	$(BINDIR)/exp1

.PHONY : run_exp
run_exp : $(BINDIR)/exp
	$(BINDIR)/exp

.PHONY : run_exp1
run_exp1 : $(BINDIR)/exp1
	$(BINDIR)/exp1

.PHONY : profile_exp
profile_exp : $(BINDIR)/exp
	$(BINDIR)/exp; gprof $(BINDIR)/exp

$(OBJDIR)/exp.o : $(OBJDIR)/orrb.o
$(OBJDIR)/exp1.o : $(OBJDIR)/orrb.o
