
TARGETS=rnd.fa rnd.k39.ctx rnd.clean.k39.ctx

CTX=../../bin/ctx63
SEQRND=../../libs/seq_file/bin/seqrnd

all: $(TARGETS)

clean:
	rm -rf $(TARGETS)

rnd.fa:
	(echo '>rnd'; $(SEQRND) 200) > rnd.fa

rnd.k39.ctx: rnd.fa
	$(CTX) build -k 39 --sample rnd --seq rnd.fa $@

rnd.clean.k39.ctx: rnd.k39.ctx
	$(CTX) clean -o $@ $<
	$(CTX) check $@

.PHONY: all clean
