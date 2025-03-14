load(file = "output/pop_gen/relatedness_all.RData")

load(file = "output/load/pheno_loads_lifetime.RData")
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf
method="gerp45"
load <- load_per_region %>% filter(loadtype == method)

# changes names
## metadata on filenames and ids
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

#merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))
gen <- left_join(gen, idnames, by = c("FID1" = "V1"))
gen <- left_join(gen, idnames, by = c("FID2" = "V1"))

#take only related
related <- subset(gen, criteria != "Unrelated")
ids_related <- data.frame(id = c(related$id.x, related$id.y)) %>% unique()
#related %>% select(c(id.x, id.y, criteria)) %>% View()

# add lms
gen <- left_join(gen, pheno_wide[,c("id", "LMS_min", "born")], by = c("id.x" = "id"))
gen <- left_join(gen, pheno_wide[,c("id", "LMS_min", "born")], by = c("id.y" = "id"))
gen <- left_join(gen, load[,c("id", "total_load")], by = c("id.x" = "id"))
gen <- left_join(gen, load[,c("id", "total_load")], by = c("id.y" = "id"))
names(gen)
gen <- gen %>% mutate(LMS_dif_abs = abs(LMS_min.x - LMS_min.y),
                      total_load_dif_abs = abs(total_load.x - total_load.y))
gen$criteria <- factor(gen$criteria, levels = c("Unrelated", "Full-sibling", "Parent-offspring", "Second-degree", "Third-degree", "Unknown"))
summary(lmerTest::lmer(LMS_dif_abs ~ criteria + (1|id.x) + (1|id.y), data = gen))
summary(lmerTest::lmer(total_load_dif_abs ~ criteria+ (1|id.x) + (1|id.y), data = gen))

parent_off <-  subset(gen, criteria == "Parent-offspring")
#choose parent
parent_off <- parent_off %>% mutate(id_parent = case_when(
  born.x > born.y ~ id.x,
  born.y > born.x ~ id.y
))

data <-  left_join(pheno_wide, load_per_region, by = "id")
`%!in%` = Negate(`%in%`)

summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), data = subset(data, id %!in% parent_off$id_parent & loadtype == "gerp45"),
                         family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), data = subset(data, id %!in% parent_off$id_parent & loadtype == "gerp45_promoter"),
                         family = "poisson", ziformula = ~1))

summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), data = subset(data, id %!in% parent_off$id_parent & loadtype == "high"),
                         family = "poisson", ziformula = ~1))

summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), data = subset(data, id %!in% parent_off$id_parent & loadtype == "high_promoter"),
                         family = "poisson", ziformula = ~1))

ids_parent_off <- data.frame(id = c(parent_off$id.x, parent_off$id.y)) %>% unique()
