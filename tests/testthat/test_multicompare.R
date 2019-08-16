context("multi.compare work")
test_that("multi.compare work", {
	data <- data.frame(A=c(0.212,0.2121,0.1213,0.2123,0.1232,0.12132, 0.8112,0.8214,0.7923,0.812322,0.8121,0.71232), 
					   B=c(0.753,0.712,0.6134,0.702,0.613,0.643,0.1213,0.1431,0.15131,0.1321,0.1632,0.1423),
					   group=c(rep("Case",6),rep("Control",6)))
	kwA <- kruskal.test(A~group, data=data)
	kwB <- kruskal.test(B~group, data=data)
	kwres <- multi.compare(fun=kruskal.test,
						   data=data,feature=c("B","A"),
						   factorNames="group")
	wxA <- wilcox.test(A~group, data=data)
	wxB <- wilcox.test(B~group, data=data)
	wxres <- multi.compare(fun=wilcox.test,data=data, 
						   feature=c("A","B"),
						   factorNames="group")
	expect_equal(kwA, kwres[[2]])
	expect_equal(kwB, kwres[[1]])
	expect_equal(wxA, wxres[[1]])
	expect_equal(wxB, wxres[[2]])
})
