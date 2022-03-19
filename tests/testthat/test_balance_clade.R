context("calculate the balance score of internal nodes")

gm.mean <- function(x, pseudonum = 0, na.rm = TRUE){
    y <- exp(mean(log(x + pseudonum), na.rm = na.rm)) - pseudonum
    return(y)
}


set.seed(123)
tr <- ape::rtree(3)
dat <- data.frame(A=c(1, 2, 3), B=c(4, 2, 1), C=c(3, 6, 5))
rownames(dat) <- tr$tip.label
mpse <- MPSE(assays=list(Abundance=dat), otutree=as.treedata(tr))

tbl1 <- mpse %>% mp_balance_clade(
           .abundance = Abundance, 
           force = T, 
           relative = F, 
           balance_fun = 'mean',
           pseudonum = 0,
           action = 'only'
         )

tbl2 <- mpse %>% mp_balance_clade(
           .abundance = Abundance, 
           force = T, 
           relative = F, 
           balance_fun = 'median',
           pseudonum = 0,
           action = 'only'
         )     

tbl3 <- mpse %>% mp_balance_clade(
           .abundance = Abundance,
           force = T,
           relative = F,
           balance_fun = 'geometric.mean',
           pseudonum = 0,
           action = 'only'
         )

binary2nodes <- tr %>% extract_binary_offspring(.node=c(4, 5))

Node1.up <- tr %>% as_tibble %>% 
        dplyr::filter(node %in% binary2nodes[[1]][[1]]) %>% 
        dplyr::pull(.data$label)

Node1.down <- tr %>% as_tibble %>%
        dplyr::filter(node %in% binary2nodes[[1]][[2]]) %>%
        dplyr::pull(.data$label)


Node2.up <- tr %>% as_tibble %>%
        dplyr::filter(node %in% binary2nodes[[2]][[1]]) %>% 
        dplyr::pull(.data$label)

Node2.down <- tr %>% as_tibble %>%
        dplyr::filter(node %in% binary2nodes[[2]][[2]]) %>% 
        dplyr::pull(.data$label)

Node1.mean <- log((dat %>% magrittr::extract(Node1.up, ) %>% apply(., 2, mean))/(dat %>% magrittr::extract(Node1.down, ) %>% apply(.,2, mean)))

Node1.median <- log((dat %>% magrittr::extract(Node1.up, ) %>% apply(., 2, median))/(dat %>% magrittr::extract(Node1.down, ) %>% apply(.,2, median)))

Node1.gm.mean <- log((dat %>% magrittr::extract(Node1.up, ) %>% apply(., 2, gm.mean))/(dat %>% magrittr::extract(Node1.down, ) %>% apply(.,2, gm.mean)))

Node2.mean <- log((dat %>% magrittr::extract(Node2.up, ) %>% apply(., 2, mean))/(dat %>% magrittr::extract(Node2.down, ) %>% apply(.,2, mean)))

Node2.median <- log((dat %>% magrittr::extract(Node2.up, ) %>% apply(., 2, median))/(dat %>% magrittr::extract(Node2.down, ) %>% apply(.,2, median)))

Node2.gm.mean <- log((dat %>% magrittr::extract(Node2.up, ) %>% apply(., 2, gm.mean))/(dat %>% magrittr::extract(Node2.down, ) %>% apply(.,2, gm.mean)))


test_that("balance score with mean",{
    expect_equal(tbl1 %>% 
           unnest(BalanceByAbundanceBySample) %>%
           filter(label == 'Node1') %>%
           pull(BalanceByAbundance, name = Sample),
           Node1.mean
    )  

    expect_equal(
           tbl1 %>%
           unnest(BalanceByAbundanceBySample) %>%	
           filter(label == 'Node2') %>% 
           pull(BalanceByAbundance, name = Sample),
           Node2.mean
    )		
  }
)


test_that("balance score with median",{
    expect_equal(tbl2 %>%
           unnest(BalanceByAbundanceBySample) %>%
           filter(label == 'Node1') %>%
           pull(BalanceByAbundance, name = Sample),
           Node1.median
    )

    expect_equal(
        tbl2 %>%
        unnest(BalanceByAbundanceBySample) %>%
        filter(label == 'Node2') %>%
        pull(BalanceByAbundance, name = Sample),
        Node2.median
    )
  }
)

test_that("balance score with geometric mean",{
    expect_equal(tbl3 %>%
                 unnest(BalanceByAbundanceBySample) %>%
                 filter(label == 'Node1') %>%
                 pull(BalanceByAbundance, name = Sample),
                 Node1.gm.mean
    )

    expect_equal(
        tbl3 %>%
        unnest(BalanceByAbundanceBySample) %>%
        filter(label == 'Node2') %>%
        pull(BalanceByAbundance, name = Sample),
        Node2.gm.mean
    )
  }
)


