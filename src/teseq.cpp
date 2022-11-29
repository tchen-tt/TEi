#include "teseq.h"

//[[Rcpp::export]]
void processbam(std::string bamfile, std::string outfq, double quantile, int length, int tsd) {
  const char *file = bamfile.data();
  const char *outfile = outfq.data();
  bam1_t *b = NULL; 
  bam_hdr_t *header = NULL;
  htsFile *fp = NULL;
  uint32_t *cigar;
  
  fp = sam_open(file, "rb");
  header = sam_hdr_read(fp);
  FILE *of = fopen(outfile, "w");
  b = bam_init1();
  
  while(sam_read1(fp, header, b) >= 0) {
    bam1_core_t *c = &b->core;
    if(c->n_cigar == 1 || (int)(c->qual) < quantile || c->flag & BAM_FDUP) continue;
    bool pass = false;
    cigar = bam_get_cigar(b);
    for(int i=0; i < c->n_cigar; ++i) {
      if (bam_cigar_op(cigar[i]) == 4 && bam_cigar_oplen(cigar[i]) > length) 
        pass = true;
    }
    //if (pass) processSequence(b, header, of, length);
    if (pass) processSequence1(b, header, of, length, tsd);
  }
  
  bam_destroy1(b);
  sam_close(fp);
  fclose(of);
}


void processSequence(bam1_t *b, bam_hdr_t *hdr, FILE *f, int length) {
  bam1_core_t *c = &b->core; 
  int32_t pos = c->pos; 
  int32_t index = 0;
  int op, op_len;
  uint32_t *cigar = bam_get_cigar(b);
  
  for (int i=0;  i < c->n_cigar; ++i) {
    op = bam_cigar_op(cigar[i]);
    op_len = bam_cigar_oplen(cigar[i]);
    if ((op != 4) || op_len < length) {
      if (bam_cigar_type(op)&2) pos += op_len;
      if (bam_cigar_type(op)&1) index += op_len;
    } 
    else {
        fprintf(f, "@%s:%d\n", hdr->target_name[c->tid], pos);
        for (int j=0; j < op_len; j++) {
          fprintf(f, "%c", int2char[bam_seqi(bam_get_seq(b), j+index)]);
        }
        fprintf(f, "\n+\n");
        for (int j=0; j < op_len; j++) {
          fprintf(f, "%c", bam_get_qual(b)[j+index] + 33);
        }
        fprintf(f, "\n");
    }
  }
}


//[[Rcpp::export]]
void processSam2bed(std::string alignmentfile, std::string outbedfile, float ratio) {
  const char *aligment = alignmentfile.data();
  const char *outbed = outbedfile.data();
  FILE *OF = fopen(outbed, "w");
  
  bam_hdr_t *header = NULL;
  bam1_t **reads = NULL, *read = NULL;
  
  htsFile *ifile = NULL;
  char *lname = NULL;
  
  int m_stack, l_stack = 0, i;
  int rv = 0;
  
  ifile = sam_open(aligment, "rb");
  header = sam_hdr_read(ifile);
  
  m_stack = 100;
  reads = (bam1_t **)malloc(sizeof(bam1_t *) * m_stack);
  
  if (reads == NULL) {
    Rcerr << "Coundn't allocate enough memery for reads array. Out of memory." << std::endl;
    rv = -2;
    free(reads);
    sam_close(ifile);
  }

  read = bam_init1();

  while(sam_read1(ifile, header, read) >= 0) {
    if(read->core.flag & BAM_FUNMAP) continue;
    if(!lname) lname = strdup(bam_get_qname(read));
    
    if(read->core.flag & BAM_FUNMAP) {
      continue;
    }
   
    if (strcmp(lname, bam_get_qname(read)) == 0) {

      if (l_stack >= m_stack) {
        m_stack *= 2;
        if ((reads = (bam1_t **)realloc(reads, sizeof(bam1_t *)* m_stack)) == NULL) {
          Rcerr << "Couldn't reallocate enough memory to hold " << 
            m_stack << "." << std::endl;
          for (i = 0; i < m_stack; i++) bam_destroy1(reads[i]);
        }
        for(i = l_stack; i<m_stack; i++) reads[i] = bam_init1();

      } 
      reads[l_stack] = bam_init1();
      bam_copy1(reads[l_stack], read);
      l_stack++;

      
    } else {
      bam2bed(reads, l_stack, OF, ratio, header);
      // bam2bed(reads, l_stack, OF, ratio, header, length);
      free(lname);
      lname = strdup(bam_get_qname(read));
      bam_copy1(reads[0], read);
      l_stack = 1;
    }
  }

  
  bam2bed(reads, l_stack, OF, ratio, header);

  for(int i=0; i<l_stack; i++) bam_destroy1(reads[i]);
  free(reads);
  sam_close(ifile);
  fclose(OF);
}

int sortFunction(const void *x, const void *y) {
  int32_t x1 = *((int32_t *)x);
  int32_t y1 = *((int32_t *)y);
  if (x1 < y1) return -1;
  if (x1 == y1) return 0;
  return 1;
}

int countNumber(int32_t *tids, int l) {
  int total = 1, i;
  int32_t last = tids[0];
  for(i = 1; i < l; i++) {
    if (last != tids[i]) {
      total++;
      last = tids[i];
    }
  }
  return total;
}

int countTable(int32_t *tids, int pos, int l) {
  double out = 1.0;
  int i;
  int32_t ltid = tids[pos];
  for(i=pos+1; i<l; i++) {
    if(tids[i] != ltid) return out;
    out += 1.0;
  }
  return out;
}

void bam2bed(bam1_t **read, int l, FILE *f, double threshold, bam_hdr_t *header) {
  int32_t tids[l], *tid_id;
  int i, elements, pos = 0, pos2, n_above = 0;
  double *counts;
  char *p;
  
  for (i = 0; i<l; i++) tids[i] = read[i]->core.tid;
  qsort(tids, l, sizeof(int32_t), sortFunction);
  elements = countNumber(tids, l);
  tid_id = (int32_t *)malloc(sizeof(int32_t) * elements);
  counts = (double *)malloc(sizeof(double) * elements);
  for (i=0; i<elements; i++) { 
    counts[i] = countTable(tids, pos, l) / ((double) l);
    tid_id[i] = tids[pos];
    pos += counts[i];
    if(counts[i] > threshold) n_above++;
  }
  if (n_above == 1) {
    for (i = 0; i < elements; i++) if(counts[i]>threshold) break;
    /*
    p = strtok(bam_get_qname(read[0]), ":");
    fprintf(f, "%s\t", p);
    p += strlen(p) + 1;
    pos2 = atoi(p);
     */
    //fprintf(f, "%i\t%i\t%s\n", (pos2-length < 0) ? 0:pos2-length, pos2+length, header->target_name[tid_id[i]]);
    //fprintf(f, "%i\t%i\t%s\t%i\t%.2f\n", (pos2-length < 0) ? 0:pos2-length, pos2+length, header->target_name[tid_id[i]], l, counts[i]);
    p = strtok(bam_get_qname(read[0]), ":");
    char_vector qname_s;
    while (p) {
      qname_s.push_back(p);
      p = strtok(NULL, ":");
    }
    fprintf(f, "%s\t", qname_s[0]);
    pos2 = atoi(qname_s[1]);
    // int length = read[0]->core.l_qseq;
    
    /*
    if (qname_s.size() == 3) {
      int read_type = atoi(qname_s[2]);
      if (read_type == 1) {
        fprintf(f, "%i\t%i\t", (pos2-length < 0) ? 0: pos2-length, pos2);
      } else {
        fprintf(f, "%i\t%i\t", pos2, pos2+length);
      }
    } else {
      fprintf(f, "%i\t%i\t", (pos2-length <0) ? 0:pos2-length, pos2+length);
    }
    fprintf(f, "%s\t", header->target_name[tid_id[i]]);
    
    if (qname_s.size() == 3) {
      fprintf(f, "%i\t", atoi(qname_s[2]));
    }
    fprintf(f, "%i\t%0.2f\n", l, counts[i]);
    */
    if (qname_s.size() == 4) {
      int read_type = atoi(qname_s[3]);
      if (read_type == 1) {
        fprintf(f, "%i\t%i\t", pos2, pos2+1);
      } 
      else {
        fprintf(f, "%i\t%i\t", pos2-1, pos2);
      }
    }
    fprintf(f, "%s\t", header->target_name[tid_id[i]]);
    if (qname_s.size() == 4) {
      fprintf(f, "%s\t", qname_s[2]);
      fprintf(f, "%i\t", atoi(qname_s[3]));
    }
    fprintf(f, "%i\t%0.2f\n", l, counts[i]);
    
  }
  free(tid_id);
  free(counts);
}


void processSequence1(bam1_t *b, bam_hdr_t *hdr, FILE *f, int length, int tsd) {
  bam1_core_t *c = &b->core; 
  int32_t pos = c->pos; 
  int32_t index = 0;
  int op, op_len;
  uint32_t *cigar = bam_get_cigar(b);
  
  for (int i=0; i < c->n_cigar; ++i) {
    op = bam_cigar_op(cigar[i]);
    op_len = bam_cigar_oplen(cigar[i]);
    if ((op != 4) || op_len < length) {
      if (bam_cigar_type(op)&2) pos += op_len;
      if (bam_cigar_type(op)&1) index += op_len;
    } 
    else {
      int type = 0;
      if (i == 0) {
        if (bam_cigar_op(cigar[i+1]) == 0) {
          type = 1; //SM
        } else {
          continue;
        }
      } 
      else if (i == c->n_cigar){
        if (bam_cigar_op(cigar[i-1]) == 0) {
          type = 2; //MS
        } else {
          continue;
        }
      } 
      else {
        if (bam_cigar_op(cigar[i-1]) == 0) {
          type = 2;
        } 
        else if (bam_cigar_op(cigar[i+1]) == 0) {
          type = 1;
        } 
        else {
          continue;
        }
      }
      
      std::string TSD;
      int32_t end;
      if (type == 1) {
        end = tsd > bam_cigar_oplen(cigar[i+1]) ? bam_cigar_oplen(cigar[i+1]) : tsd;
        for (int j = 0; j < end; ++j) {
          char n = int2char[bam_seqi(bam_get_seq(b), j+op_len+index)];
          TSD.push_back(n);
        }
      } else if (type == 2) {
        end = tsd > bam_cigar_oplen(cigar[i-1]) ? bam_cigar_oplen(cigar[i-1]) : tsd;
        for (int j = 0; j < end; ++j) {
          char n = int2char[bam_seqi(bam_get_seq(b), index-j-1)];
          TSD.push_back(n);
          //std::reverse(TSD.begin(), TSD.end());
        }
        std::reverse(TSD.begin(), TSD.end());
      }
      
      fprintf(f, "@%s:%d:%s:%i\n", hdr->target_name[c->tid], pos, TSD.data(),type);
        
      for (int j=0; j < op_len; j++) {
        fprintf(f, "%c", int2char[bam_seqi(bam_get_seq(b), j+index)]);
      }
      fprintf(f, "\n+\n");
      for (int j=0; j < op_len; j++) {
        fprintf(f, "%c", bam_get_qual(b)[j+index] + 33);
      }
      fprintf(f, "\n");
      
    }
  }
}

