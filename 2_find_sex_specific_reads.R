# parse a fasta file of unmapped reads, extract the ones that are male specific

# libraries
library("ape")
library("dplyr")

# fasta file

unmapped.fa <- read.dna(file = "unmapped_reads.fa", format = "fasta", as.character = TRUE, as.matrix = FALSE)

# id sex file

id.sex <- read.table("id_sex.txt")
names(id.sex) <- c("id", "sex")

# locus names
locus.names <- names(unmapped.fa)

#split locus labels into glorious dataframe

fasta1 <- unmapped.fa[1]

options(error = NULL)

fasta_to_row <- function(label){
  
  fasta <- unmapped.fa[label]
  id <- unlist(strsplit(label, split = " "))[2] %>% gsub("\\[|\\]", "", .)
  c.locus <- unlist(strsplit(label, split = "_"))[2] %>% as.numeric()
  allele <- unlist(strsplit(label, split = " "))[1] %>% strsplit(.,split = "_") %>% unlist %>% .[8] %>% as.numeric()
  sequence <- fasta[[1]] %>% toupper %>% paste0(collapse = "") 
  out.df <- data.frame(id, c.locus, allele, sequence)
  return(out.df)

}

unmapped.df <- lapply(names(unmapped.fa), fasta_to_row)
unmapped.df <- bind_rows(unmapped.df)
write.table(unmapped.df, file = "unmapped.df.txt", quote = FALSE, row.names = FALSE)
unmapped.df <- read.table(file = "unmapped.df.txt", header = TRUE, stringsAsFactors = FALSE)

# add in sex info and arrange
unmapped.df <- left_join(id.sex, unmapped.df)

unmapped.df <- unmapped.df %>%
  arrange(c.locus, sex, id)

# find that loci that occur in males
loci.males <- unmapped.df %>%
  filter(sex == "M") %>%
  filter(!is.na(c.locus)) %>%
  select(c.locus) %>%
  unique %>%
  .[,1]

loci.males %>% head 


# find that loci that occur in females
loci.females <- unmapped.df %>%
  filter(sex == "F") %>%
  filter(!is.na(c.locus)) %>%
  select(c.locus) %>%
  unique %>%
  .[,1]


# find that loci that *only* occur in males
loci.male.only <- loci.males [!(loci.males %in% loci.females)]

male.loci.df <- unmapped.df %>%
  filter(c.locus %in% loci.male.only)

# read in species data

species.dat <- read.table(file = "whtcmn_faststructure.txt", header = TRUE, stringsAsFactors = FALSE)
species.dat <- species.dat %>%
  filter(k.value.run == 2) %>%
  select(id, membership)

species.dat$id <- species.dat$id %>% paste0("_2014")

# join in cluster data for individuals who have it 
male.loci.df <- left_join(male.loci.df, species.dat)
male.loci.df$pop <- substr(male.loci.df$id, 1, 2)

# use population assignments to approximate species for individuals missing cluster data
fake_memberships <- function (pop){
  if (pop %in% c("WR", "PQ", "CP", "SK", "BR", "AL", "MR")){
    return(1)
  }else{
    return(2)
  }
  
}

male.loci.df$membership[is.na(male.loci.df$membership)] <- lapply(male.loci.df$pop[is.na(male.loci.df$membership)], fake_memberships) %>% unlist
  
# c.loci spefici to whtstbk?
# result: nope

loci.white <- male.loci.df %>%
  filter(membership == 2) %>%
  filter(!is.na(c.locus)) %>%
  select(c.locus) %>%
  unique %>%
  .[,1]

loci.common <- male.loci.df %>%
  filter(membership == 1) %>%
  filter(!is.na(c.locus)) %>%
  select(c.locus) %>%
  unique %>%
  .[,1]

loci.white.only <- loci.white [!(loci.white %in% loci.common)]

male.loci.df %>%
  filter(c.locus %in% loci.white.only)

# sequences specific to the white stickleback?

